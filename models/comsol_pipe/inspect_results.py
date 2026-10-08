"""Extract conserved fluxes, field profiles and a preview from a solved MPH."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import numpy as np
from .build_model import COMSOL, HERE
from models.config import COMSOL_TEMPLATE, DATA_DIR
from .helium_properties import boundary_diagnostics


def inspect(model, destination):
    destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=True)
    j=model.java
    co=j.component('comp1')
    for tag, typ, selection in [('maxgas','Maximum','geom1_gas_dom'),('mingas','Minimum','geom1_gas_dom')]:
        if tag not in list(co.cpl().tags()):
            op=co.cpl().create(tag,typ); op.selection().named(selection)
    v=co.variable('diagnostics')
    for key, expression in {
        'divU':'d(u,r)+if(r>0[m],u/r,d(u,r))+d(w,z)',
        'tauzz':'2*muHe(T)*d(w,z)-2/3*muHe(T)*divU',
        'taurz':'muHe(T)*(d(u,z)+d(w,r))',
        'Fenergy':'rhoHe*w*(cpHe*T+0.5*(u^2+w^2))-u*taurz-w*tauzz',
    }.items():
        v.set(key,expression)
    # Refresh postprocessing definitions added to an already computed model.
    j.sol('sol1').updateSolution()
    exprs={
        'mass_in_kg_s':'comp1.intin(comp1.rhoHe*w)',
        'mass_out_kg_s':'comp1.intout(comp1.rhoHe*w)',
        'cold_outward_conduction_W':'comp1.intcold(comp1.ht.ndflux)',
        'warm_outward_conduction_W':'comp1.intwarm(comp1.ht.ndflux)',
        'steel_cold_outward_conduction_W':'comp1.intsc(comp1.ht.ndflux)',
        'steel_warm_outward_conduction_W':'comp1.intsw(comp1.ht.ndflux)',
        'gas_advective_stress_energy_in_W':'comp1.intin(comp1.Fenergy)',
        'gas_advective_stress_energy_out_W':'comp1.intout(comp1.Fenergy)',
        'comsol_energy_balance_W':'comp1.ht.energyBalance',
        'comsol_heat_balance_W':'comp1.ht.heatBalance',
        'max_Kn':'comp1.maxgas(comp1.Kn)',
        'max_Mach':'comp1.maxgas(comp1.Mach)',
        'max_pointwise_Re':'comp1.maxgas(comp1.Re_local)',
        'min_pressure_Pa':'comp1.mingas(p)',
        'max_pressure_Pa':'comp1.maxgas(p)',
        'min_temperature_K':'comp1.mingas(T)',
        'max_temperature_K':'comp1.maxgas(T)',
        'max_speed_m_s':'comp1.maxgas(comp1.spf.U)',
        'gas_inlet_area_m2':'comp1.intin(1)',
        'gas_outlet_area_m2':'comp1.intout(1)',
        'mean_inlet_pressure_Pa':'comp1.intin(p)/comp1.intin(1)',
        'mean_outlet_pressure_Pa':'comp1.intout(p)/comp1.intout(1)',
        'kSS_at_1K_W_mK':'kSS(1[K])',
        'kSS_at_293K_W_mK':'kSS(293.15[K])',
        'kHe_at_1K_W_mK':'kHe(1[K])',
    }
    results={'model_name':model.name(),'continuum_baseline_only':True,'values':{},'evaluation_errors':{}}
    # Couplings made through the Java API use unweighted meridional measure.
    # Apply the cylindrical measure once, including in average denominators.
    for key,expr in list(exprs.items()):
        for opname in ('intin','intout','intcold','intwarm','intsc','intsw'):
            expr=expr.replace(f'comp1.{opname}(',f'comp1.{opname}(2*pi*r*')
        exprs[key]=expr
    for key,expr in exprs.items():
        try:
            val=float(model.evaluate(expr))
            results['values'][key]=val
            print(key,val,flush=True)
        except Exception as error:
            results['evaluation_errors'][key]=str(error)
            print('EVALUATION ERROR',key,str(error)[:500],flush=True)
    vals=results['values']
    if 'mass_in_kg_s' in vals and 'mass_out_kg_s' in vals:
        vals['relative_mass_imbalance']=abs(vals['mass_out_kg_s']-vals['mass_in_kg_s'])/max(abs(vals['mass_in_kg_s']),abs(vals['mass_out_kg_s']))
    flux_keys=['cold_outward_conduction_W','warm_outward_conduction_W','gas_advective_stress_energy_out_W','gas_advective_stress_energy_in_W']
    if all(k in vals for k in flux_keys):
        a,b,c,d=[vals[k] for k in flux_keys]
        vals['independent_total_energy_residual_W']=a+b+c-d
        vals['relative_total_energy_residual']=abs(a+b+c-d)/max(abs(a),abs(b),abs(c),abs(d))
    # Store a reviewable derived-values table inside the MPH, not only in JSON.
    if 'summary' in list(j.result().numerical().tags()):
        j.result().numerical().remove('summary')
    if 'tblSummary' not in list(j.result().table().tags()):
        j.result().table().create('tblSummary','Table')
    j.result().table('tblSummary').label('Verified fluxes and continuum diagnostics')
    summary=j.result().numerical().create('summary','EvalGlobal')
    summary.label('Mass, heat, pressure and validity diagnostics')
    summary.set('expr',list(exprs.values()))
    summary.set('descr',list(exprs.keys()))
    summary.set('table','tblSummary')
    summary.setResult()
    (destination/'results.json').write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding='utf-8')

    # A common z sampling includes near-end gradients as well as the full tube.
    z=np.unique(np.r_[np.linspace(0,.288,401),np.geomspace(1e-7,.02,120),.288-np.geomspace(1e-7,.02,120)])
    out=[]
    for radial,label,expressions in [(0.,'axis',['T','p','w','Kn','Mach']),(.003175-.000381,'wall',['T'])]:
        if 'profile_'+label in list(j.result().numerical().tags()):
            j.result().numerical().remove('profile_'+label)
        ev=j.result().numerical().create('profile_'+label,'Interp')
        ev.set('expr',expressions)
        ev.set('coord',np.array([np.full_like(z,radial),z]).tolist())
        values=np.asarray(ev.getData(),dtype=float).squeeze()
        if values.ndim==1: values=values[None,:]
        np.savetxt(destination/f'{label}_profile.csv',np.column_stack([z,values.T]),delimiter=',',header=','.join(['z_m']+expressions),comments='')
        out.append(values)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10,6.7),layout='constrained')
    gas,wall=out
    axes[0,0].plot(z*1000,wall[0],label='Steel wall',lw=2)
    axes[0,0].plot(z*1000,gas[0],label='Gas centerline',ls='--')
    axes[0,0].set_ylabel('Temperature (K)'); axes[0,0].legend()
    axes[0,1].plot(z*1000,gas[1]); axes[0,1].set_ylabel('Centerline absolute pressure (Pa)')
    axes[1,0].plot(z*1000,gas[2]); axes[1,0].set_ylabel('Centerline axial velocity (m/s)')
    axes[1,1].semilogy(z*1000,gas[3]); axes[1,1].axhline(.01,color='gray',ls=':',label='Kn = 0.01')
    axes[1,1].set_ylabel('Centerline Knudsen number'); axes[1,1].legend()
    for ax in axes.flat:
        ax.set_xlabel('Distance from cold inlet (mm)'); ax.grid(alpha=.2)
    fig.suptitle('He4 flow in insulated SS304 pipe — stationary continuum baseline\nL = 288 mm, bore = 4.826 mm; rarefied warm region limits physical accuracy',fontsize=12)
    fig.savefig(destination/'profiles.png',dpi=180)
    plt.close(fig)
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model',type=Path,nargs='?',default=COMSOL_TEMPLATE)
    parser.add_argument('--output',type=Path,default=DATA_DIR/'pipe_inspection')
    args=parser.parse_args()
    os.environ['PATH']=str(COMSOL.parent)+os.pathsep+os.environ['PATH']
    import mph
    client=mph.Client(cores=4,version='6.1')
    model=client.load(str(args.model.resolve()))
    if model.java.sol('sol1').isEmpty():
        raise RuntimeError('Model has no computed solution. Run build_model without --build-only.')
    inspect(model,args.output)
    model.save(str(args.model.resolve()))


if __name__=='__main__': main()
