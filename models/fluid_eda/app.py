"""Local schematic editor and asynchronous, cancellable simulation API."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import json
from pathlib import Path
import threading
import time
from urllib.parse import urlsplit
import uuid
from flask import Flask, jsonify, render_template, request, send_from_directory
from .examples import continuous_pump_example, DEFAULTS
from .network import validate_network, NetworkError
from .simulation import run_simulation
from .error_reporting import describe_error
from .comsol_backend import COMSOL_EXECUTABLE,TEMPLATE,ComsolCancelled

from models.config import DATA_DIR


def create_app(*,results_dir=None,runner=run_simulation):
    app=Flask(__name__)
    app.config.update(MAX_CONTENT_LENGTH=1024*1024,JSON_AS_ASCII=False)
    app.json.ensure_ascii=False
    root=Path(results_dir or DATA_DIR/'jobs');root.mkdir(parents=True,exist_ok=True)
    jobs={};lock=threading.RLock();executor=ThreadPoolExecutor(max_workers=1,thread_name_prefix='fluid-eda')
    app.extensions['eda']={'jobs':jobs,'executor':executor,'results_dir':root}

    @app.before_request
    def local_requests_only():
        if request.host.split(':')[0] not in ('127.0.0.1','localhost','[::1]'):
            return jsonify(error={'message':'Only local requests are accepted','code':'local_only'}),403
        origin=request.headers.get('Origin')
        if origin and (urlsplit(origin).netloc!=request.host or urlsplit(origin).scheme!='http'):
            return jsonify(error={'message':'Cross-origin requests are not accepted','code':'cross_origin'}),403

    @app.after_request
    def headers(response):
        response.headers['X-Content-Type-Options']='nosniff'
        response.headers['Cache-Control']='no-store' if request.path.startswith('/api/') else 'no-cache'
        response.headers['Referrer-Policy']='same-origin'
        return response

    @app.get('/')
    def index():return render_template('index.html')

    @app.get('/api/example')
    def example():return jsonify(continuous_pump_example())

    @app.get('/api/catalog')
    def catalog():return jsonify(DEFAULTS)

    @app.get('/api/status')
    def status():
        return jsonify(status='ready',comsol_available=COMSOL_EXECUTABLE.is_file() and TEMPLATE.is_file(),
                       comsol={'available':COMSOL_EXECUTABLE.is_file() and TEMPLATE.is_file(),'version':'6.1'},
                       engines=['reduced','comsol','two_phase'],version='0.1.0')

    def snapshot(job):
        return {k:v for k,v in job.items() if k not in ('cancel_event','future','project')}

    def execute(job):
        job_id=job['job_id'];folder=root/job_id
        def progress(fraction,message):
            with lock:job.update(progress=max(job['progress'],min(1,float(fraction))),message=str(message))
        try:
            if job['cancel_event'].is_set():raise ComsolCancelled('Calculation cancelled')
            with lock:job.update(status='running',message='Starting circuit solver')
            folder.mkdir(parents=True,exist_ok=True)
            (folder/'project.json').write_text(json.dumps(job['project'],ensure_ascii=False,indent=2),encoding='utf-8')
            result=runner(job['project'],progress=progress,cancel_event=job['cancel_event'],output_dir=folder)
            if job['cancel_event'].is_set():raise ComsolCancelled('Calculation cancelled')
            # Reject nonfinite results at the boundary rather than ship invalid JSON.
            payload=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)
            (folder/'results.json').write_text(payload,encoding='utf-8')
            with lock:job.update(status='complete',progress=1.,message='Calculation complete',result=result)
        except ComsolCancelled:
            with lock:job.update(status='cancelled',message='Calculation cancelled')
        except Exception as error:
            detail=describe_error(error,job['project'],stage='Circuit solve')
            with lock:job.update(status='error',message=detail['message'],error=detail)
            app.logger.exception('Simulation failed: %s',job_id)
        finally:
            with lock:
                job['elapsed_seconds']=time.monotonic()-job.pop('_started',time.monotonic())
                data=snapshot(job)
            folder.mkdir(parents=True,exist_ok=True)
            (folder/'job.json').write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

    @app.post('/api/simulate')
    def simulate():
        project=request.get_json(silent=True)
        try:
            validate_network(project)
            if project.get('settings',{}).get('engine','comsol') not in ('comsol','reduced','two_phase'):
                raise NetworkError('Unknown pipe-solver mode')
            json.dumps(project,allow_nan=False)
        except (ValueError,TypeError) as error:
            detail=describe_error(error,project,code='invalid_json')
            return jsonify(error=detail),400
        with lock:
            active=sum(j['status'] in ('queued','running') for j in jobs.values())
            if active>=3:return jsonify(error={'message':'Calculations are already queued, Wait for or cancel the current job','code':'queue_full'}),429
            jid=uuid.uuid4().hex
            job={'job_id':jid,'status':'queued','progress':0.,'message':'Waiting for solver',
                 'created_at':datetime.now(timezone.utc).isoformat(),'_started':time.monotonic(),
                 'result':None,'error':None,'project':project,'cancel_event':threading.Event()}
            jobs[jid]=job;job['future']=executor.submit(execute,job)
        return jsonify(job_id=jid),202

    @app.get('/api/jobs/<job_id>')
    def get_job(job_id):
        with lock:
            if job_id not in jobs:return jsonify(error={'message':'Simulation job not found','code':'not_found'}),404
            data=snapshot(jobs[job_id]);data.pop('_started',None)
            return jsonify(data)

    @app.post('/api/jobs/<job_id>/cancel')
    def cancel(job_id):
        with lock:
            job=jobs.get(job_id)
            if not job:return jsonify(error={'message':'Simulation job not found'}),404
            if job['status'] in ('running','queued'):
                job['cancel_event'].set();job['message']='Cancelling calculation…'
                if job['future'].cancel():job.update(status='cancelled',message='Calculation cancelled')
            return jsonify(snapshot(job))

    @app.get('/api/jobs/<job_id>/files/<path:filename>')
    def files(job_id,filename):
        if job_id not in jobs or jobs[job_id]['status']!='complete':return jsonify(error={'message':'Results are not yet available'}),404
        return send_from_directory(root/job_id,filename,as_attachment=True)

    @app.errorhandler(413)
    def too_large(_):return jsonify(error={'message':'Circuit file is too large (limit 1 MB)','code':'too_large'}),413
    return app


def main():
    parser=argparse.ArgumentParser(description='Helium Studio local fluid EDA')
    parser.add_argument('--port',type=int,default=5004);args=parser.parse_args()
    print(f'Helium Studio: http://127.0.0.1:{args.port}',flush=True)
    create_app().run(host='127.0.0.1',port=args.port,debug=False,threaded=True,use_reloader=False)

if __name__=='__main__':main()
