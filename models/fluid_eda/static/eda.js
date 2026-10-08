/* Helium Studio: standalone schematic editor; all stored model values use SI. */
(() => {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const svg = $('schematic');
  const NS = 'http://www.w3.org/2000/svg';
  const STORE = 'helium-studio-project-v1';
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const clone = (v) => JSON.parse(JSON.stringify(v));
  const forwardBackflowDefaults = {forward_speed_m3_s:.240,compression_ratio:1e8,forward_speed_curve:[],compression_curve:[],max_inlet_pressure_Pa:null,max_outlet_pressure_Pa:null,max_throughput_Pa_m3_s:null};
  const icons = {
    pipe: '<path d="M1 8h7m24 0h7M1 20h7m24 0h7M8 5h24v18H8Z M11 5v18m18-18v18"/>',
    pump: '<circle cx="20" cy="14" r="11"/><path d="M1 14h8m22 0h8M16 8l10 6-10 6Z"/>',
    source: '<circle cx="20" cy="14" r="11"/><path d="M1 14h8m22 0h8M13 17q3-4 7 0t7 0M17 8v5m6-5v5"/>',
    jt: '<path d="M1 14h10m18 0h10M11 5v18l18-18v18L11 5Z"/>',
    junction: '<path d="M3 14h34M20 2v24"/><circle cx="20" cy="14" r="3" fill="currentColor"/>',
  };
  const types = {
    pipe: {name:'Stainless-steel pipe',en:'SS 304 PIPE',prefix:'PIPE',width:120,height:38,params:{length_m:.288,outer_diameter_m:.00635,wall_thickness_m:.000762,T_a_K:293.15,T_b_K:293.15,wall_boundary:'adiabatic',wall_profile:[]}},
    pump: {name:'Vacuum pump',en:'VACUUM PUMP',prefix:'PUMP',width:94,height:62,params:{mode:'constant_speed',speed_m3_s:.004,expression:'Pin*S/(R*Tin)',outlet_temperature_K:293.15,...forwardBackflowDefaults}},
    source: {name:'Pressure source / Pot',en:'PRESSURE SOURCE',prefix:'POT',width:94,height:66,params:{mode:'constant',pressure_Pa:15.57,temperature_K:1,expression:'15.57',heat_load_W:0}},
    jt: {name:'J–T Restriction',en:'J–T RESTRICTION',prefix:'JT',width:90,height:38,params:{law:'capillary',material:'CuNi',diameter_m:.0001,length_m:.1,impedance_mode:'measured',impedance_m_inv3:3.2e16,discharge_coefficient:.7}},
    junction: {name:'Junction',en:'JUNCTION',prefix:'J',width:16,height:16,params:{temperature_K:293.15}},
  };
  let project = freshProject();
  const state = {selected:null,tab:'params',tool:'select',viewport:{x:0,y:0,w:1080,h:650},initialWidth:1080,undo:[],redo:[],pending:null,pointer:null,space:false,canvasFocused:false,result:null,resultRevision:null,revision:0,job:null,jobRevision:null,pollTimer:null,toastTimer:null,saveTimer:null,profile:'temperature',color:'temperature',status:null,validationErrors:[],errorMessage:'',error:null};

  function freshProject() { return {version:1,name:'Untitled circuit',settings:{engine:'comsol'},components:[],wires:[]}; }
  function toast(message, error=false, duration=4000) { clearTimeout(state.toastTimer); $('toast').textContent=message; $('toast').classList.toggle('error',error); $('toast').hidden=false; state.toastTimer=setTimeout(()=>{$('toast').hidden=true;},duration); }
  function commit() { state.undo.push(JSON.stringify(project)); if(state.undo.length>60) state.undo.shift(); state.redo=[]; }
  function clearError() {state.error=null;state.validationErrors=[];state.errorMessage='';}
  function changed({render=true}={}) {const hadError=!!state.error;state.revision++;clearError();autosave();renderResults();if(render)renderAll();else{renderScene();if(hadError)renderInspector();}updateToolbar();}
  function autosave() { $('save-state').textContent='Saving…'; clearTimeout(state.saveTimer); state.saveTimer=setTimeout(()=>{try{localStorage.setItem(STORE,JSON.stringify(project));$('save-state').textContent='Saved locally';}catch(_){$('save-state').textContent='Export to save';}},250); }
  function resetResult() {state.result=null;state.resultRevision=null;renderResults();}
  function updateToolbar() { $('undo-btn').disabled=!state.undo.length;$('redo-btn').disabled=!state.redo.length;$('delete-btn').disabled=!state.selected;$('rotate-btn').disabled=!selectedComponent()||selectedComponent().type==='junction';$('flip-btn').disabled=$('rotate-btn').disabled;$('component-count').textContent=`${project.components.length} components`;$('canvas-state').textContent=state.pending?'Connecting':state.job?'Calculating':state.selected?'Selected':'Ready'; }
  function selectedComponent() {return project.components.find(c=>c.id===state.selected);}
  function undo() {if(!state.undo.length)return;state.redo.push(JSON.stringify(project));project=JSON.parse(state.undo.pop());state.selected=null;state.pending=null;changed();syncControls();}
  function redo() {if(!state.redo.length)return;state.undo.push(JSON.stringify(project));project=JSON.parse(state.redo.pop());state.selected=null;state.pending=null;changed();syncControls();}
  function syncControls() {$('project-name').value=project.name || 'Untitled circuit';$('engine-select').value=project.settings?.engine || 'comsol';}
  function uid(prefix='c') {return `${prefix}_${Date.now().toString(36)}_${Math.random().toString(36).slice(2,7)}`;}
  function snap(v) {return Math.round(v/20)*20;}
  function addComponent(type, position) {
    if(!types[type])return;
    commit();const t=types[type];let n=1;while(project.components.some(c=>c.label===`${t.prefix}-${String(n).padStart(2,'0')}`))n++;
    const pos=position || {x:state.viewport.x+state.viewport.w/2,y:state.viewport.y+state.viewport.h/2};
    let x=snap(pos.x),y=snap(pos.y);if(!position){while(project.components.some(c=>Math.abs(c.x-x)<90&&Math.abs(c.y-y)<85)){x+=40;y+=80;}}
    const c={id:uid(),type,label:`${t.prefix}-${String(n).padStart(2,'0')}`,x,y,rotation:0,flipX:false,params:clone(t.params)};
    project.components.push(c);state.selected=c.id;state.tab='params';changed();svg.focus();state.canvasFocused=true;
  }
  function removeSelected() { if(!state.selected)return;commit();const id=state.selected;project.components=project.components.filter(c=>c.id!==id);project.wires=project.wires.filter(w=>w.id!==id && w.from.component!==id && w.to.component!==id);state.selected=null;state.pending=null;changed(); }
  function rotateSelected() {const c=selectedComponent();if(!c || c.type==='junction')return;commit();c.rotation=((c.rotation||0)+90)%360;changed();}
  function flipSelected() {const c=selectedComponent();if(!c || c.type==='junction')return;commit();c.flipX=!c.flipX;changed();}
  function duplicateSelected() {const c=selectedComponent();if(!c)return;commit();const copy=clone(c);copy.id=uid();copy.label=`${copy.label} Copy`;copy.x+=40;copy.y+=80;project.components.push(copy);state.selected=copy.id;changed();}
  function select(id, tab) {state.selected=id;if(tab)state.tab=tab;renderScene();renderInspector();renderResultsTable();updateToolbar();}

  function worldPoint(event) {const pt=svg.createSVGPoint();pt.x=event.clientX;pt.y=event.clientY;return pt.matrixTransform(svg.getScreenCTM().inverse());}
  // Mirror in screen coordinates after rotation; port identities never change.
  function localPoint(c,x,y=0) {const theta=(c.rotation||0)*Math.PI/180;return {x:(c.flipX?-1:1)*(x*Math.cos(theta)-y*Math.sin(theta)),y:x*Math.sin(theta)+y*Math.cos(theta)};}
  function portPosition(component,port) {const t=types[component.type];const p=localPoint(component,component.type==='junction'?0:port==='a'?-t.width/2:t.width/2);return {x:component.x+p.x,y:component.y+p.y};}
  function endpoint(ep) {const c=project.components.find(c=>c.id===ep.component);return c?portPosition(c,ep.port):null;}
  function wirePath(a,b,from,to) {
    const component=from && project.components.find(c=>c.id===from.component);
    const vertical=component && (component.rotation||0)%180!==0;
    const other=to && project.components.find(c=>c.id===to.component);
    if(component && other && component.type!=='junction' && other.type!=='junction') {
      const dir=(c,p)=>{const point=localPoint(c,p==='a'?-1:1);return {x:Math.round(point.x),y:Math.round(point.y)};};
      const d1=dir(component,from.port),d2=dir(other,to.port),s={x:a.x+d1.x*24,y:a.y+d1.y*24},e={x:b.x+d2.x*24,y:b.y+d2.y*24};
      let points=[a,s];
      if(d1.y===0 && d2.y===0){
        if((d1.x!==d2.x&&d1.x*(e.x-s.x)<0)||(d1.x===d2.x&&Math.abs(s.y-e.y)<48)){const y=Math.abs(s.y-e.y)>=48?(s.y+e.y)/2:Math.min(s.y,e.y)-80;points.push({x:s.x,y},{x:e.x,y});}
        else{const x=d1.x===d2.x?(d1.x>0?Math.max(s.x,e.x)+24:Math.min(s.x,e.x)-24):(s.x+e.x)/2;points.push({x,y:s.y},{x,y:e.y});}
      }
      else if(d1.x===0 && d2.x===0){
        if((d1.y!==d2.y&&d1.y*(e.y-s.y)<0)||(d1.y===d2.y&&Math.abs(s.x-e.x)<48)){const x=Math.abs(s.x-e.x)>=48?(s.x+e.x)/2:Math.min(s.x,e.x)-100;points.push({x,y:s.y},{x,y:e.y});}
        else{const y=d1.y===d2.y?(d1.y>0?Math.max(s.y,e.y)+24:Math.min(s.y,e.y)-24):(s.y+e.y)/2;points.push({x:s.x,y},{x:e.x,y});}
      }
      else points.push(d1.x===0?{x:s.x,y:e.y}:{x:e.x,y:s.y});
      points.push(e,b);return points.map((p,i)=>`${i?'L':'M'}${p.x},${p.y}`).join('');
    }
    if(Math.abs(a.x-b.x)<2 || Math.abs(a.y-b.y)<2)return `M${a.x},${a.y}L${b.x},${b.y}`;
    if(vertical){const y=(a.y+b.y)/2;return `M${a.x},${a.y}L${a.x},${y}L${b.x},${y}L${b.x},${b.y}`;}
    const x=(a.x+b.x)/2;return `M${a.x},${a.y}L${x},${a.y}L${x},${b.y}L${b.x},${b.y}`;
  }
  function connectPort(id,port) {
    if(!state.pending){state.pending={component:id,port};$('connection-hint').hidden=false;$('canvas-wrap').classList.add('connecting');renderScene();updateToolbar();return;}
    const from=state.pending,to={component:id,port};
    if(from.component===to.component){if(from.port===to.port)cancelConnection();else toast('Select a port on another component.');return;}
    const key=(v)=>`${v.component}:${v.port}`;
    if(project.wires.some(w=>(key(w.from)===key(from)&&key(w.to)===key(to))||(key(w.from)===key(to)&&key(w.to)===key(from)))){toast('These ports are already connected.');cancelConnection();return;}
    commit();project.wires.push({id:uid('w'),from:clone(from),to});cancelConnection();changed();
  }
  function cancelConnection() {state.pending=null;$('pending-wire').setAttribute('hidden','');$('connection-hint').hidden=true;$('canvas-wrap').classList.remove('connecting');renderScene();updateToolbar();}
  function componentSubtitle(c) {
    const p=c.params;
    if(c.type==='pipe')return `${fmt(p.length_m*1000,4)} mm · Ø${fmt(p.outer_diameter_m*1000,4)}`;
    if(c.type==='pump')return p.mode==='forward_backflow'?'Forward–backflow':p.mode==='expression'?'Custom flow curve':`${fmt(p.speed_m3_s*1000,4)} L/s`;
    if(c.type==='source')return p.mode==='saturation'?`${fmt(p.temperature_K,4)} K · Saturated vapor`:p.mode==='expression'?'P = f(T)':`${fmt(p.pressure_Pa,4)} Pa · ${fmt(p.temperature_K,4)} K`;
    if(c.type==='jt')return `${p.law==='capillary'?`${p.material||''} capillary`:'Orifice'} · Ø${fmt(p.diameter_m*1000,4)} mm${p.law==='capillary'?` × ${fmt(p.length_m*1000,4)} mm`:''}`;
    return 'Mixing node';
  }
  function symbol(c) {
    const w=types[c.type].width;
    if(c.type==='pipe')return `<path class="symbol" d="M-60,-11H-49M49,-11H60M-60,11H-49M49,11H60"/><rect class="body" x="-49" y="-17" width="98" height="34" rx="3"/><path class="symbol" d="M-42,-17V17M42,-17V17"/><path class="symbol" opacity=".45" d="M-22,0H22m-5,-4 5,4-5,4"/>`;
    if(c.type==='pump')return `<path class="symbol" d="M-${w/2},0H-29M29,0H${w/2}"/><circle class="body" r="29"/><path class="symbol" d="M-10,-15L17,0L-10,15Z"/>`;
    if(c.type==='source')return `<path class="symbol" d="M-${w/2},0H-29M29,0H${w/2}"/><circle class="body" r="29"/><path class="symbol" d="M-20,7Q-13,0-6,7T8,7T22,7M-10,-16V-6M0,-19V-9M10,-16V-6"/>`;
    if(c.type==='jt')return `<path class="symbol" d="M-45,0H-19M19,0H45"/><path class="body" d="M-19,-16V16L19,-16V16Z"/><path class="symbol" d="M0,-16V-26m-6,0H6"/>`;
    return '<circle class="body" r="6"/>';
  }
  function getComponentResult(id) {const c=state.result?.components;return Array.isArray(c)?c.find(v=>v.id===id):c?.[id];}
  function colorDomain() {const values=project.components.map(c=>getComponentResult(c.id)).filter(Boolean).flatMap(r=>state.color==='pressure'?[r.pressure_in,r.pressure_out]:[r.temperature_in,r.temperature_out]).filter(Number.isFinite);return values.length?[Math.min(...values),Math.max(...values)]:null;}
  function colorScale(value,domain) {if(!domain||!Number.isFinite(value))return null;const t=domain[1]===domain[0]?.5:Math.max(0,Math.min(1,(value-domain[0])/(domain[1]-domain[0])));const a=[64,135,201],b=[216,144,86];return `rgb(${a.map((n,i)=>Math.round(n+(b[i]-n)*t)).join(',')})`;}
  function renderScene() {
    $('connection-hint').hidden=!state.pending;$('canvas-wrap').classList.toggle('connecting',!!state.pending);if(!state.pending)$('pending-wire').setAttribute('hidden','');
    const domain=state.color!=='none'?colorDomain():null;
    $('wire-layer').innerHTML=project.wires.map(w=>{const a=endpoint(w.from),b=endpoint(w.to);if(!a||!b)return '';const path=wirePath(a,b,w.from,w.to);return `<g class="wire ${state.selected===w.id?'selected':''}" data-wire="${esc(w.id)}"><path class="wire-hit" d="${path}"/><path class="wire-path" d="${path}"/></g>`;}).join('');
    $('component-layer').innerHTML=project.components.map(c=>{
      const t=types[c.type];if(!t)return '';const r=getComponentResult(c.id);const resultValue=r?(state.color==='pressure'?(r.pressure_in+r.pressure_out)/2:(r.temperature_in+r.temperature_out)/2):null;const color=colorScale(resultValue,domain);const angle=c.rotation||0;
      const ports=c.type==='junction'?['a']:['a','b'];
      const invalid=state.validationErrors.includes(c.id);
      const portMarkup=ports.map(p=>{const x=c.type==='junction'?0:p==='a'?-t.width/2:t.width/2,pos=localPoint(c,x),label=localPoint(c,x+(p==='a'?-8:8),-9);const pending=state.pending?.component===c.id&&state.pending?.port===p;const failed=invalid&&state.error?.context?.port===p;return `<g data-port="${p}" data-component="${esc(c.id)}" role="button" aria-label="${esc(c.label)} Port ${p}${failed?' — simulation error':''}"><circle class="port-hit" cx="${pos.x}" cy="${pos.y}" r="12"/><circle class="port ${pending?'connecting':''} ${failed?'port-error':''}" cx="${pos.x}" cy="${pos.y}" r="${failed?6:4}"/>${c.type!=='junction'?`<text class="port-label" x="${label.x}" y="${label.y}" text-anchor="middle">${p.toUpperCase()}</text>`:''}</g>`;}).join('');
      const isVertical=angle%180!==0,labelY=isVertical?-7:t.height/2+23,labelX=isVertical?t.height/2+20:0,anchor=isVertical?'start':'middle';
      return `<g class="component ${state.selected===c.id?'selected':''} ${invalid?'component-error':''}" data-id="${esc(c.id)}" data-testid="component-${esc(c.type)}" transform="translate(${c.x},${c.y})"><title>${esc(c.label)} · ${esc(t.name)}</title><g class="component-symbol" transform="scale(${c.flipX?-1:1},1) rotate(${angle})"><rect class="selection-box" x="${-t.width/2-15}" y="${-t.height/2-12}" width="${t.width+30}" height="${t.height+24}"/>${symbol(c)}${color?`<path d="M${-Math.min(24,t.width/3)},${t.height/2+6}H${Math.min(24,t.width/3)}" stroke="${color}" stroke-width="3" stroke-linecap="round"/>`:''}</g>${portMarkup}${invalid?`<g class="error-badge" transform="translate(${isVertical?0:t.width/2+10},${isVertical?-t.width/2-20:-t.height/2-17})"><title>Simulation error: ${esc(state.errorMessage)}</title><circle r="9"/><text text-anchor="middle" y="4">!</text></g>`:''}<text class="component-label" x="${labelX}" y="${labelY}" style="text-anchor:${anchor}">${esc(c.label)}</text><text class="component-value" x="${labelX}" y="${labelY+16}" style="text-anchor:${anchor}">${esc(componentSubtitle(c))}</text>${r&&Number.isFinite(r.mass_flow)?`<text class="result-badge" x="${labelX}" y="${labelY+32}" style="text-anchor:${anchor}">${fmt(r.mass_flow*1e9,3)} µg/s</text>`:''}</g>`;
    }).join('');
    $('canvas-empty').hidden=project.components.length>0;
    $('color-legend').hidden=!domain;
    if(domain){$('legend-min').textContent=fmt(domain[0],3);$('legend-max').textContent=fmt(domain[1],3);$('legend-label').textContent=state.color==='pressure'?'Pressure / Pa':'Temperature / K';}
    updateViewport();
  }
  function updateViewport() {const v=state.viewport;svg.setAttribute('viewBox',`${v.x} ${v.y} ${v.w} ${v.h}`);$('zoom-reset').textContent=`${Math.round(state.initialWidth/v.w*100)}%`;}
  function zoom(factor,center) {const v=state.viewport;const newW=Math.max(280,Math.min(9000,v.w*factor));factor=newW/v.w;const p=center||{x:v.x+v.w/2,y:v.y+v.h/2};v.x=p.x-(p.x-v.x)*factor;v.y=p.y-(p.y-v.y)*factor;v.w*=factor;v.h*=factor;updateViewport();}
  function fit() {const box=svg.getBoundingClientRect();if(!box.width||!box.height)return;const ratio=box.width/box.height;let minX=0,minY=0,maxX=1080,maxY=650;if(project.components.length){minX=Math.min(...project.components.map(c=>c.x))-120;maxX=Math.max(...project.components.map(c=>c.x))+120;minY=Math.min(...project.components.map(c=>c.y))-105;maxY=Math.max(...project.components.map(c=>c.y))+145;}
    let w=Math.max(maxX-minX,450),h=Math.max(maxY-minY,320);if(w/h>ratio)h=w/ratio;else w=h*ratio;state.viewport={x:(minX+maxX-w)/2,y:(minY+maxY-h)/2,w,h};state.initialWidth=w;updateViewport();}
  function setTool(tool) {state.tool=tool;$('select-tool').classList.toggle('active',tool==='select');$('pan-tool').classList.toggle('active',tool==='pan');$('canvas-wrap').classList.toggle('pan',tool==='pan');}

  function field(label,key,value,unit='',scale=1,extra={}) {return `<label class="field"><span class="field-label">${label}${extra.required?'<span class="pill-unit">Boundary conditions</span>':''}</span><span class="input-unit"><input data-param="${key}" data-scale="${scale}" aria-label="${label}" type="number" value="${Number.isFinite(value)?Number((value*scale).toPrecision(12)):''}" step="any" ${extra.min!==undefined?`min="${extra.min}"`:''} ${extra.max!==undefined?`max="${extra.max}"`:''}><span>${unit}</span></span></label>`;}
  function selectField(label,key,value,options) {return `<label class="field">${label}<select data-param="${key}" aria-label="${label}">${options.map(([v,l])=>`<option value="${v}" ${v===value?'selected':''}>${l}</option>`).join('')}</select></label>`;}
  function help(title,body) {return `<details class="help field-help"><summary><span class="help-icon">?</span>${title}</summary><p>${body}</p></details>`;}
  function section(title,body) {return `<div class="form-section"><div class="form-section-title">${title}</div>${body}</div>`;}
  function renderInspector() {
    $('params-tab').classList.toggle('active',state.tab==='params');$('outputs-tab').classList.toggle('active',state.tab==='outputs');$('params-tab').setAttribute('aria-selected',state.tab==='params');$('outputs-tab').setAttribute('aria-selected',state.tab==='outputs');
    const target=$('inspector-content'),c=selectedComponent();
    if(!c){if(state.selected && project.wires.some(w=>w.id===state.selected)){target.innerHTML=`<div class="device-heading"><span class="device-icon">⌁</span><div><h2>Wire</h2><span class="eyebrow">IDEAL CONNECTION</span></div></div><p class="muted" style="font-size:11px;line-height:1.9">Connected ports share one node pressure. Wires have no flow resistance.</p><button class="button secondary" id="inspector-delete">Delete connection</button>`;$('inspector-delete').onclick=removeSelected;return;}
      target.innerHTML='<div class="inspector-empty"><svg viewBox="0 0 48 48"><rect x="8" y="9" width="32" height="30" rx="5"/><path d="M14 18h20M14 24h12M14 30h17"/><circle cx="33" cy="30" r="2"/></svg><h3>Edit every detail</h3><p>Select a component to inspect its geometry,<br>boundary conditions and simulation results.</p></div>';return;
    }
    const t=types[c.type],p=c.params;
    const heading=`<div class="device-heading"><span class="device-icon"><svg viewBox="0 0 40 28">${icons[c.type]}</svg></span><div><h2>${esc(c.label)}</h2><span class="eyebrow">${t.en}</span></div></div>${state.validationErrors.includes(c.id)?`<div class="validation-error"><b>Simulation error${state.error?.context?.port?` · Port ${esc(state.error.context.port.toUpperCase())}`:''}</b><p>${esc(state.errorMessage)}</p>${state.error?.hint?`<p>${esc(state.error.hint)}</p>`:''}</div>`:''}`;
    if(state.tab==='outputs'){target.innerHTML=heading+resultInspector(c);bindProfile();return;}
    let content=`<label class="field">Component name<input id="component-name" aria-label="Component name" value="${esc(c.label)}" maxlength="50"></label>`;
    if(c.type==='pipe'){
      content+=section('Geometry · Stainless steel 304',field('Pipe length','length_m',p.length_m,'mm',1000,{min:.001})+`<div class="field-grid">${field('Outer diameter','outer_diameter_m',p.outer_diameter_m,'mm',1000,{min:.001})}${field('Wall thickness','wall_thickness_m',p.wall_thickness_m,'mm',1000,{min:.00001})}</div>`);
      content+=section('Thermal boundary',selectField('Outer-wall boundary','wall_boundary',p.wall_boundary||'adiabatic',[['adiabatic','Adiabatic outer wall'],['prescribed','Piecewise outer-wall temperature']])+(p.wall_boundary==='prescribed'?'':`<div class="field-grid">${field('A end temperature','T_a_K',p.T_a_K,'K',1,{min:.01})}${field('B end temperature','T_b_K',p.T_b_K,'K',1,{min:.01})}</div>`)+help('Pipe boundaries and heat conduction',project.settings.engine==='two_phase'?'Specify wall temperatures. Fluid energy is advected through finite gas-wall exchange; fluid outlet temperature is computed. Axial steel conduction is retained in the adiabatic suction pipe. Conduction between components is excluded.':p.wall_boundary==='prescribed'?'Outer-wall temperature follows the piecewise-linear profile below. Prescribed cooling stages set the local wall temperature. End temperatures follow the first and last segments. Gas and steel-wall heat transfer are solved separately; the gas temperature need not equal the prescribed wall temperature.':'Gas and wall temperatures are fixed at each end. The outer wall is adiabatic. The model includes gas-wall heat exchange and axial conduction. The model retains conduction within this pipe, but end conduction is not transferred to neighboring components or added to the pot energy balance.'));
      if(p.wall_boundary==='prescribed') content+=wallProfileEditor(p);
      if(project.settings.engine==='two_phase')content+=field('Heat-transfer multiplier','heat_transfer_multiplier',p.heat_transfer_multiplier??1,'',1,{min:.01,max:1000})+help('Two-phase heat transfer','Finite exchange determines fluid temperature and quality. The starting estimate is Nu = 3.66 with homogeneous mixture transport; multiplier 1 is a model assumption, not measured condensation data. Use a calibrated multiplier or sensitivity sweep. Wall endpoints prescribe steel temperature; fluid outlet temperature is solved.');
      content+=help('Pipe solver mode',project.settings.engine==='two_phase'?'Experimental 1D continuous-pot mode: one pot, suction pipe, pump, precooler and capillary connected A to B. Public real-helium properties and a common circulation flow are solved together.':project.settings.engine==='comsol'?'COMSOL solves each parameterized 2D axisymmetric pipe using the network end pressures.':'Quick preview uses a cross-section-averaged 1D model. Select COMSOL in the left panel to compute full pipe solutions.');
    } else if(c.type==='pump'){
      content+=selectField('Pump characteristics','mode',p.mode,[['constant_speed','Constant pumping speed'],['expression','Custom mass-flow function'],['forward_backflow','Forward–backflow']]);
      if(p.mode==='forward_backflow')content+=forwardBackflowEditor(p);
      else {
        content+=field('Pumping speed','speed_m3_s',p.speed_m3_s,'L/s',1000,{min:0});
        if(p.mode==='expression') content+=`<label class="field">ṁ = f(Pin, Pout, Tin)<textarea data-param="expression" aria-label="Pump mass-flow expression" spellcheck="false">${esc(p.expression)}</textarea><span class="field-note">Output units kg/s</span></label>`;
        content+=help('Pump flow definition',p.mode==='expression'?'Variables: Pin and Pout (Pa), Tin (K), S (m³/s), and R (helium specific gas constant). Return the A → B mass flow in kg/s.':'Pumping speed S is the inlet volumetric flow rate. Mass flow is ṁ = Pin · S / (R · Tin), positive from A to B.');
      }
      content+=section('Discharge boundary',field('Discharge temperature','outlet_temperature_K',p.outlet_temperature_K,'K',1,{min:.01})+help('Temperature boundary','Specify the port B temperature after compression and heat exchange with the pump. Pump work and cooling power are not calculated separately. Forward–backflow uses the port A temperature to convert both directional throughputs to mass flow; unequal port temperatures require calibration.'));
    } else if(c.type==='source'){
      content+=selectField('Pressure constraint','mode',p.mode,[['constant','Constant pressure'],['saturation','He-4 Saturated vapor pot'],['expression','Temperature function P(T)']]);
      content+=field('Node temperature','temperature_K',p.temperature_K,'K',1,{min:.01});
      if(p.mode==='constant')content+=field('Fixed pressure','pressure_Pa',p.pressure_Pa,'Pa',1,{min:.000001});
      if(p.mode==='expression')content+=`<label class="field">P = f(T)<textarea data-param="expression" aria-label="Pressure-temperature expression" spellcheck="false">${esc(p.expression)}</textarea><span class="field-note">T in K; output in Pa</span></label>`;
      if(p.mode==='saturation')content+=field('Applied heat load','heat_load_W',p.heat_load_W??0,'mW',1000,{min:0})+help('Applied heat load','Sample and other heat entering the pot, excluding intercomponent conduction. Default zero is a no-load trial. Match this load to net cooling when trialing the equilibrium temperature.');
      content+=help('Pressure source and pot',p.mode==='saturation'?'Node pressure follows the He-4 saturated vapor-pressure correlation at the specified temperature. The specified temperature is a trial pot temperature. Incoming helium enthalpy is retained across the JT restriction; saturated vapor leaves the pot. Results compare net cooling with the applied heat load. Solid conduction between components is excluded.':'Ports A and B share one fixed-pressure node. The source supplies or absorbs net mass as required by the network. Outgoing gas has the specified temperature.');
    } else if(c.type==='jt'){
      content+=selectField('Flow-resistance model','law',p.law,[['orifice','Compressible orifice · includes choking'],['capillary','Laminar capillary']]);
      content+=field(p.law==='capillary'?'Capillary inner diameter':'Orifice diameter','diameter_m',p.diameter_m,'mm',1000,{min:.000001});
      if(p.law==='capillary'){
        content+=field('Capillary length','length_m',p.length_m,'mm',1000,{min:.000001});
        content+=`<label class="field">Material<input data-param="material" aria-label="Restriction material" value="${esc(p.material||'')}" maxlength="80"></label>`;
        content+=selectField('Capillary resistance','impedance_mode',p.impedance_mode||'geometry',[['measured','Measured impedance'],['geometry','Calculate from geometry']]);
        if(p.impedance_mode==='measured')content+=field('Measured impedance z','impedance_m_inv3',p.impedance_m_inv3,'m⁻³',1,{min:1e-30});
        content+=help('Impedance definition','z = ΔP / (μ V̇), in m⁻³. Geometry gives z = 128 L / (π d⁴). Measured mode uses the entered z instead: changing dimensions does not rescale it. The gas model uses ṁ = (Pa² − Pb²) / (2 μ R T z). CuNi is a material label; restriction-wall heat conduction is not modeled.');
      }
      else content+=field('Discharge coefficient Cd','discharge_coefficient',p.discharge_coefficient,'',1,{min:.001,max:1});
      if(project.settings.engine==='two_phase')content+=field('Wall roughness','roughness_m',p.roughness_m??0,'micrometre',1e6,{min:0})+help('Capillary friction','HEM uses Churchill friction and includes acceleration and choking. Default zero roughness assumes a smooth wall. Measured z specifies effective hydraulic length through z = 128 Leff / (pi d^4); it does not independently calibrate two-phase friction.');
      content+=help('Flow resistance and throttling temperature',project.settings.engine==='two_phase'?'The capillary conserves total enthalpy including kinetic energy. Phase change affects density, friction and critical flow. The final pot flash includes dissipation of exit-jet kinetic energy. Normal-liquid transport is never extended into He II.':p.law==='capillary'?'The capillary uses the compressible laminar pressure-squared flow law, neglecting environmental heat exchange. Ideal helium has no isenthalpic temperature drop. Real J–T cooling and liquefaction require real-fluid properties.':'The orifice uses a compressible gas-flow equation with choking below the critical pressure ratio. Ideal helium has no isenthalpic temperature drop. Real J–T cooling and liquefaction require real-fluid properties.');
    } else if(c.type==='junction'){
      content+=field('Node temperature','temperature_K',p.temperature_K,'K',1,{min:.01})+help('Junction','Join multiple wires at one zero-resistance node. Its specified node temperature must match adjacent pipe-end states. Component-to-component conduction is excluded; gas flow energy remains in the balance.');
    }
    content+=`<div class="inspector-actions">${c.type!=='junction'?'<button id="inspector-flip" class="button secondary" title="Flip horizontally (X)">Flip ↔</button>':''}<button id="duplicate-btn" class="button secondary">Duplicate</button><button id="inspector-delete" class="button subtle danger">Delete</button></div>`;
    target.innerHTML=heading+content;
    $('component-name').addEventListener('change',event=>{const value=event.target.value.trim();if(!value){event.target.value=c.label;return;}commit();c.label=value;changed({render:false});});
    target.querySelectorAll('[data-param]').forEach(input=>input.addEventListener('change',()=>{
      const key=input.dataset.param;let value=input.value;
      if(input.type==='number'){value=Number(value)/(Number(input.dataset.scale)||1);if(input.dataset.optional && input.value.trim()==='')value=null;else if(!Number.isFinite(value)||input.value.trim()===''||!input.checkValidity()){toast('Enter a valid parameter value.',true);input.value=c.params[key]===null?'':String(c.params[key]*(Number(input.dataset.scale)||1));return;}}
      commit();c.params[key]=value;
      if(key==='wall_boundary' && value==='prescribed' && !c.params.wall_profile?.length){c.params.wall_profile=[{length_m:c.params.length_m,T_start_K:c.params.T_a_K,T_end_K:c.params.T_b_K}];}
      changed({render:false});if(input.tagName==='SELECT')renderInspector();
    }));
    bindWallProfile(c);
    bindPumpCurves(c);
    if($('inspector-flip'))$('inspector-flip').onclick=flipSelected;$('duplicate-btn').onclick=duplicateSelected;$('inspector-delete').onclick=removeSelected;
  }
  function forwardBackflowEditor(p) {
    const curve=(label,key,scale,columns)=>`<label class="field">${label}<textarea data-pump-curve="${key}" data-scale="${scale}" aria-label="${label}" spellcheck="false" placeholder="${columns}">${esc((p[key]||[]).map(([pressure,value])=>`${pressure}, ${Number((value*scale).toPrecision(12))}`).join('\n'))}</textarea></label>`;
    const optional=(label,key,unit)=>field(label,key,p[key],unit,1,{min:1e-15}).replace('data-param=', 'data-optional="true" data-param=');
    return section('Forward and backward flow',
      field('Nominal forward speed','forward_speed_m3_s',p.forward_speed_m3_s,'L/s',1000,{min:1e-12})+
      field('Zero-flow compression K₀','compression_ratio',p.compression_ratio,'',1,{min:1})+
      help('Equation and HiPace data','ṁ = Sf(Pin) · [Pin − Pout / K₀(Pout)] / (R · Tin). Positive flow is A → B; negative flow is B → A. There is no zero-flow clamp. Both terms use the inlet-temperature reference. The <a href="https://www.pfeiffervacuum.com/ca/en/product/hipace-300-m.html" target="_blank" rel="noopener">HiPace 300 M helium specifications</a> give 240 L/s and K₀ = 10⁸. These nominal defaults are an approximation, not a measured two-pressure map. The TM700 manual lists 215 L/s and K₀ &gt; 10⁸; confirm your hardware revision.')+
      '<button class="button subtle small" id="pump-hipace-defaults">Use HiPace nominal data</button>')+
      `<details class="help field-help pump-advanced"><summary><span class="help-icon">?</span>Pressure-dependent curves</summary><div class="pump-advanced-content">`+
      curve('Forward-speed table · Pa, L/s','forward_speed_curve',1000,'Inlet pressure (Pa), forward speed (L/s)')+
      curve('Compression table · Pa, K₀','compression_curve',1,'Outlet pressure (Pa), compression ratio')+
      '<p>Optional: one pair per line, at least two rows, with positive pressures in ascending order. Leave blank to use the nominal constant. Tables override the corresponding nominal value. Interpolation is log–log; values outside a table hold the nearest endpoint and produce a result warning. Measured net speed approximates forward speed only when reference backflow is negligible.</p></div></details>'+
      `<details class="help field-help pump-advanced"><summary><span class="help-icon">?</span>Operating limits</summary><div class="pump-advanced-content">`+
      optional('Maximum inlet pressure','max_inlet_pressure_Pa','Pa')+
      optional('Maximum outlet pressure','max_outlet_pressure_Pa','Pa')+
      optional('Maximum |net throughput|','max_throughput_Pa_m3_s','Pa·m³/s')+
      '<p>Leave unknown limits blank. Enter verified helium limits for your pump, speed and cooling conditions. The converged operating point is checked against these limits. Nitrogen forevacuum ratings are not helium limits. Throughput uses the inlet-temperature reference.</p></div></details>';
  }
  function validPumpCurve(rows,key,scale=1) {
    return Array.isArray(rows)&&(rows.length===0||rows.length>=2&&rows.length<=128)&&rows.every((row,i)=>
      Array.isArray(row)&&row.length===2&&row.every(Number.isFinite)&&row[0]>=1e-9&&row[0]<=1e9&&
      (i===0||row[0]>rows[i-1][0])&&row[1]/scale>=(key==='compression_curve'?1:1e-15)&&row[1]/scale<=(key==='compression_curve'?1e20:1e6));
  }
  function bindPumpCurves(c) {
    if(c.type!=='pump'||c.params.mode!=='forward_backflow')return;
    $('inspector-content').querySelectorAll('[data-pump-curve]').forEach(input=>{
      input.oninput=()=>input.setCustomValidity('');
      input.onchange=()=>{
        const key=input.dataset.pumpCurve,scale=Number(input.dataset.scale);
        const rows=input.value.trim()?input.value.trim().split(/\r?\n/).filter(line=>line.trim()).map(line=>line.trim().split(/[,\s]+/).map(Number)):[];
        if(!validPumpCurve(rows,key,scale)){input.setCustomValidity('Use 2–128 pressure/value pairs, with strictly increasing pressures in 1e-9–1e9 Pa. Speed must be positive; compression must be 1–1e20.');input.reportValidity();return;}
        input.setCustomValidity('');commit();c.params[key]=rows.map(([pressure,value])=>[pressure,value/scale]);changed({render:false});
      };
    });
    $('pump-hipace-defaults').onclick=()=>{commit();Object.assign(c.params,{forward_speed_m3_s:.240,compression_ratio:1e8,forward_speed_curve:[],compression_curve:[]});changed();toast('HiPace nominal data applied. Entered operating limits were retained.');};
  }
  function fmt(n,d=4) {if(n===null||n===undefined||!Number.isFinite(Number(n)))return '—';n=Number(n);if(n===0)return '0';if(Math.abs(n)>=1e6 || Math.abs(n)<.0001)return n.toExponential(Math.max(1,d-1));return Number(n.toPrecision(d)).toLocaleString('en-US',{maximumFractionDigits:Math.max(0,d-Math.floor(Math.log10(Math.abs(n)))-1)});}
  function wallProfileEditor(p) {
    const rows=p.wall_profile||[];
    return section('Outer-wall segments',`<div class="wall-profile-header"><span>Length / mm</span><span>Start / K</span><span>End / K</span><span></span></div><div class="wall-profile-rows">${rows.map((row,index)=>`<div class="wall-profile-row"><input type="number" step="any" min="0.000001" value="${Number((row.length_m*1000).toPrecision(10))}" data-wall-row="${index}" data-wall-key="length_m" aria-label="Segment ${index+1} length in mm"><input type="number" step="any" min="0.01" value="${row.T_start_K}" data-wall-row="${index}" data-wall-key="T_start_K" aria-label="Segment ${index+1} start temperature in K"><input type="number" step="any" min="0.01" value="${row.T_end_K}" data-wall-row="${index}" data-wall-key="T_end_K" aria-label="Segment ${index+1} end temperature in K"><button data-wall-remove="${index}" class="wall-remove" title="Delete segment ${index+1}" aria-label="Delete segment ${index+1}" ${rows.length<2?'disabled':''}>×</button></div>`).join('')}</div><button class="button subtle small wall-add" id="wall-add">＋ Add segment</button><div class="wall-profile-total">Total <b>${fmt(rows.reduce((s,r)=>s+r.length_m,0)*1000,6)} mm</b><span>Pipe length updates automatically</span></div>`+help('Wall temperature profile','Wall temperature varies linearly within each segment. Equal start and end temperatures define an isothermal segment. Segment edits update the pipe length and end temperatures. These temperatures prescribe the local precooling boundary.'));
  }
  function bindWallProfile(c) {
    if(c.type!=='pipe'||c.params.wall_boundary!=='prescribed')return;
    for(const key of ['length_m','T_a_K','T_b_K']){const input=$('inspector-content').querySelector(`[data-param="${key}"]`);if(input){input.readOnly=true;input.title='Updated automatically from the wall profile';}}
    function sync(){const rows=c.params.wall_profile;c.params.length_m=rows.reduce((sum,row)=>sum+row.length_m,0);c.params.T_a_K=rows[0].T_start_K;c.params.T_b_K=rows[rows.length-1].T_end_K;changed();}
    $('inspector-content').querySelectorAll('[data-wall-row]').forEach(input=>{input.onchange=()=>{let value=Number(input.value);if(!Number.isFinite(value)||!input.value.trim()||!input.checkValidity()){toast('Segment lengths and temperatures must be positive.',true);renderInspector();return;}const key=input.dataset.wallKey;if(key==='length_m')value/=1000;commit();c.params.wall_profile[Number(input.dataset.wallRow)][key]=value;sync();};});
    $('inspector-content').querySelectorAll('[data-wall-remove]').forEach(button=>{button.onclick=()=>{if(c.params.wall_profile.length<2)return;commit();c.params.wall_profile.splice(Number(button.dataset.wallRemove),1);sync();};});
    if($('wall-add')){$('wall-add').disabled=c.params.wall_profile.length>=32;$('wall-add').onclick=()=>{if(c.params.wall_profile.length>=32)return;const last=c.params.wall_profile[c.params.wall_profile.length-1];commit();c.params.wall_profile.push({length_m:.1,T_start_K:last?.T_end_K||293.15,T_end_K:last?.T_end_K||293.15});sync();};}
  }
  function output(label,value,unit,wide=false){return `<div class="output-value ${wide?'wide':''}"><span>${label}</span><b>${fmt(value,5)}</b><small>${unit}</small></div>`;}
  function resultInspector(c) {
    const r=getComponentResult(c.id);if(!r)return '<div class="inspector-empty"><h3>Awaiting simulation results</h3><p>Run a simulation to view<br>temperature, pressure and axial profiles.</p></div>';
    let html='<div class="output-grid">'+output('Inlet pressure',r.pressure_in??r.pressure,'Pa')+output(state.result?.engine==='two_phase'&&c.type==='jt'?'Pot back pressure':'Outlet pressure',r.pressure_out??r.pressure,'Pa')+output('Inlet temperature',r.temperature_in??r.temperature,'K')+output(state.result?.engine==='two_phase'&&c.type==='jt'?'Pot-entry equilibrium T':'Outlet temperature',r.temperature_out??r.temperature,'K')+output(c.type==='source'?'Net mass supply':'Mass flow',Number.isFinite(r.mass_flow)?r.mass_flow*1e9:null,'µg/s',true)+'</div>';
    if(c.type==='source'&&c.params.mode==='saturation'){
      const e=r.pot_energy;
      if(e?.available){
        html+=section('Pot energy balance','<div class="output-grid">'+
          output('Net cooling power',e.net_cooling_power_W*1000,'mW',true)+
          output('Applied heat load',e.applied_heat_load_W*1000,'mW')+
          output('Heat imbalance',e.energy_residual_W*1000,'mW')+'</div>');
        const streams=e.incoming_streams.map(v=>`${esc(project.components.find(p=>p.id===v.component_id)?.label||v.component_id)}: ${esc(v.upstream_phase)} feed; h = ${fmt(v.enthalpy_J_kg/1000,5)} kJ/kg`+
          (Number.isFinite(v.state_at_pot_pressure.vapor_quality)?`; vapor after throttling ${fmt(v.state_at_pot_pressure.vapor_quality*100,4)}%`:''));
        html+=help('Pot balance and signs',`Evaporation term: ${fmt(e.evaporation_power_W*1000,5)} mW<br>Incoming enthalpy load: ${fmt(e.incoming_enthalpy_load_W*1000,5)} mW<br>`+
          `Inflow: ${fmt(e.mass_in_kg_s*1e9,5)} microgram/s; outflow: ${fmt(e.mass_out_kg_s*1e9,5)} microgram/s<br>`+streams.join('<br>')+
          '<br><br>Net cooling = outgoing vapor energy minus incoming feed energy, referenced to pot liquid. Positive cooling can remove a load; negative cooling heats the pot. Heat imbalance = applied load minus net cooling: positive warms, negative cools. Equilibrium requires this value and the mass imbalance to approach zero. Incoming JT enthalpy is retained. Conduction between components is excluded.'+
          '<br><br>Pot saturation properties use DB98; normal-fluid feed enthalpy uses CoolProp with an aligned reference. '+(e.model?.startsWith('HEM_')?'The coupled HEM mode solves circulation and phase change together.':'Existing pipe/restriction flow laws remain gas approximations.')+
          (e.notes?.length?'<br><br>'+e.notes.map(esc).join('<br>'):''));
      } else {
        html+=help('Pot energy unavailable',esc(e?.message||'Run the simulation again to calculate the pot enthalpy balance.'));
      }
    } else if(Number.isFinite(r.fluid_cooling_power_W)){
      html+=section('Fluid energy','<div class="output-grid">'+output('Net flow cooling',r.fluid_cooling_power_W*1000,'mW',true)+'</div>');
    }
    if(r.pump_mode==='forward_backflow'){
      const d=r.diagnostics||{},direction=d.flow_direction==='backflow'?'Backflow · B → A':d.flow_direction==='balanced'?'Balanced · zero net flow':'Forward · A → B';
      html+=section('Pump flow balance',`<div class="pump-direction">${direction}</div><div class="output-grid">`+
        output('Forward contribution',d.forward_mass_flow_kg_s*1e9,'µg/s')+output('Backflow contribution',d.backflow_mass_flow_kg_s*1e9,'µg/s')+
        output('Effective net speed',r.pumping_speed_m3_s*1000,'L/s')+output('Compression K₀',d.compression_ratio,'')+
        output('Zero-flow inlet pressure',d.zero_flow_inlet_pressure_Pa,'Pa')+output('Net throughput',d.net_throughput_Pa_m3_s,'Pa·m³/s')+'</div>');
    }
    if(state.result?.engine==='two_phase'&&r.profile?.length){
      const last=r.profile[r.profile.length-1];
      html+=section(c.type==='jt'?'Capillary exit state':'Fluid outlet state','<div class="output-grid">'+`<div class="output-value wide"><span>Phase</span><b>${esc(String(last.phase||'gas').replaceAll('_',' '))}</b></div>`+(Number.isFinite(last.vapor_quality)?output('Vapor mass fraction',last.vapor_quality*100,'%'):'' )+'</div>');
      if(c.type==='jt')html+=help('Exit and pot states',`Choked: ${r.diagnostics?.choked?'yes':'no'}. Capillary exit pressure: ${fmt(last.pressure,6)} Pa. The profile ends at the capillary exit; downstream mixing into the pot is a separate control volume. Pot-pressure state: ${esc(r.flash_at_pot?.phase||'unavailable')}${Number.isFinite(r.flash_at_pot?.vapor_quality)?`; vapor mass fraction ${fmt(r.flash_at_pot.vapor_quality*100,5)}%`:''}.`);
    }
    if(r.profile?.length>1){
      const options=[['temperature','Fluid temperature / K'],['vapor_quality','Vapor mass fraction'],['void_fraction','Vapor volume fraction'],['enthalpy_J_kg','Specific enthalpy / J/kg'],['wall_temperature','Wall temperature / K'],['outer_wall_temperature','Outer-wall temperature / K'],['pressure','Pressure / Pa'],['velocity','Velocity / m·s⁻¹']].filter(([key])=>r.profile.some(p=>Number.isFinite(p[key])));
      if(!options.some(([key])=>key===state.profile))state.profile=options[0]?.[0]||'temperature';
      html+=`<div class="chart-heading"><b>Axial profiles</b><select id="profile-select" aria-label="Profile quantity">${options.map(([key,label])=>`<option value="${key}" ${state.profile===key?'selected':''}>${label}</option>`).join('')}</select></div><svg id="profile-chart" class="profile-chart" viewBox="0 0 260 190" role="img" aria-label="Component axial profiles"></svg><div class="chart-note" id="profile-note"></div>`;
    }
    const ds={...r,...r.diagnostics};ds.max_Kn??=ds.max_knudsen??ds.max_kn;ds.max_Re??=ds.max_reynolds??ds.max_re;
    const diagnosticFields=[['max_Kn','Maximum Kn','',1],['max_Re','Maximum Re','',1]];
    const items=diagnosticFields.filter(([key])=>Number.isFinite(ds[key]));if(items.length)html+=section('Additional results',items.map(([key,label,unit,scale])=>`<div class="field-label" style="font-size:11px;margin:12px 0;color:#8a9bb2"><span>${label}</span><b style="font-weight:500;color:#5c7597">${fmt(ds[key]*scale,5)} ${unit}</b></div>`).join(''));
    if(r.warnings?.length)html+=help('Component model notes',r.warnings.map(esc).join('<br>'));
    return html;
  }
  function bindProfile(){if($('profile-select')){$('profile-select').onchange=event=>{state.profile=event.target.value;drawProfile();};drawProfile();}}
  function drawProfile(){const chart=$('profile-chart'),c=selectedComponent();if(!chart||!c)return;const data=getComponentResult(c.id)?.profile?.filter(p=>Number.isFinite(p.x)&&Number.isFinite(p[state.profile]));if(!data?.length){chart.innerHTML='<text x="130" y="90" text-anchor="middle">No data available for this quantity</text>';return;}
    const hasWall=state.profile==='temperature'&&data.some(p=>Number.isFinite(p.wall_temperature));
    if($('profile-note'))$('profile-note').textContent=`Position x / mm, A → B${hasWall?' · Solid: fluid, Dashed: wall':state.profile==='wall_temperature'?' · Steel-wall mid-thickness':state.profile==='outer_wall_temperature'?' · Outer steel surface':c.type==='pipe'?' · cross-section average':''}`;
    const ys=data.map(p=>p[state.profile]);if(hasWall)ys.push(...data.filter(p=>Number.isFinite(p.wall_temperature)).map(p=>p.wall_temperature));
    const x0=Math.min(...data.map(p=>p.x)),x1=Math.max(...data.map(p=>p.x));let y0=Math.min(...ys),y1=Math.max(...ys);let spread=y1-y0;if(spread===0)spread=Math.abs(y1)*.1||1;y0-=spread*.08;y1+=spread*.08;if(['vapor_quality','void_fraction'].includes(state.profile)){y0=0;y1=1;}
    const left=46,right=251,top=15,bottom=158;const X=x=>left+(x-x0)/(x1-x0||1)*(right-left),Y=y=>bottom-(y-y0)/(y1-y0)*(bottom-top);
    const axisDigits=Math.min(9,Math.max(3,Math.ceil(Math.log10(Math.max(Math.abs(y0),Math.abs(y1),1)))-Math.floor(Math.log10(spread))+2));
    let out='';for(let i=0;i<4;i++){const y=y0+(y1-y0)*i/3,yy=Y(y);out+=`<line class="grid" x1="${left}" y1="${yy}" x2="${right}" y2="${yy}"/><text x="${left-8}" y="${yy+3}" text-anchor="end">${fmt(y,axisDigits)}</text>`;}
    for(let i=0;i<3;i++){const x=x0+(x1-x0)*i/2;out+=`<text x="${X(x)}" y="${bottom+17}" text-anchor="middle">${fmt(x*1000,4)}</text>`;}
    const curve=data.map((p,i)=>`${i?'L':'M'}${X(p.x).toFixed(2)},${Y(p[state.profile]).toFixed(2)}`).join('');out+=`<path class="area" d="${curve}L${X(data[data.length-1].x)},${bottom}L${X(data[0].x)},${bottom}Z"/><path class="curve" d="${curve}"/>`;
    if(hasWall){const wall=data.filter(p=>Number.isFinite(p.wall_temperature)).map((p,i)=>`${i?'L':'M'}${X(p.x).toFixed(2)},${Y(p.wall_temperature).toFixed(2)}`).join('');out+=`<path d="${wall}" fill="none" stroke="#d39664" stroke-width="1.5" stroke-dasharray="4 3"/>`;}
    chart.innerHTML=out;
  }
  function renderResults() {
    const result=state.result;const running=!!state.job;renderError();
    $('result-empty').hidden=!!result||running||!!state.error;$('progress-block').hidden=!running;$('result-table-wrap').hidden=!result||running;$('result-summary').hidden=!result||running;$('download-results').disabled=!result;$('warnings').hidden=!result?.warnings?.length||running;
    if(!result&&!running&&!state.error){$('result-status').className='status-pill neutral';$('result-status').textContent='Awaiting calculation';}
    if(result&&!running){const stale=state.resultRevision!==state.revision;const warn=(result.warnings?.length||0)>0;$('result-status').textContent=stale?'Results are out of date':warn?'Complete · Model warnings':'Complete';$('result-status').className=`status-pill ${stale||warn?'warning':'success'}`;const s=result.summary||{};
      const engine=s.engine||result.engine||project.settings.engine;
      let summary=`<span>Solver mode<b>${engine==='comsol'?'COMSOL coupled simulation':engine==='reduced'?'1D quick preview':engine==='two_phase'?'Two-phase continuous pot':esc(engine)}</b></span>`;
      if(Number.isFinite(s.mass_balance_relative))summary+=`<span>Mass-balance residual<b>${fmt(s.mass_balance_relative,2)}</b></span>`;
      if(Number.isFinite(s.max_Kn))summary+=`<span>Maximum Kn<b>${fmt(s.max_Kn,3)}</b></span>`;
      $('result-summary').innerHTML=summary;
      $('warnings').innerHTML=result.warnings?.length?`<details><summary>${result.warnings.length} model warnings · Click to view</summary><ul>${result.warnings.map(w=>`<li>${esc(typeof w==='string'?w:w.message||JSON.stringify(w))}</li>`).join('')}</ul></details>`:'';renderResultsTable();
    }
    if(state.error){$('result-status').textContent='Simulation failed';$('result-status').className='status-pill error';}
  }
  function renderResultsTable(){if(!state.result)return;$('result-table').innerHTML=project.components.map(c=>{const r=getComponentResult(c.id);if(!r)return '';return `<tr data-result-id="${esc(c.id)}" class="${state.selected===c.id?'selected':''}"><td>${esc(c.label)}</td><td>${fmt(r.pressure_in??r.pressure,5)}</td><td>${fmt(r.pressure_out??r.pressure,5)}</td><td>${fmt(r.temperature_in??r.temperature,5)}</td><td>${fmt(r.temperature_out??r.temperature,5)}</td><td>${fmt(Number.isFinite(r.mass_flow)?r.mass_flow*1e9:null,5)}</td></tr>`;}).join('');}
  function renderAll(){renderScene();renderInspector();updateToolbar();}

  async function api(url,options={}) {const response=await fetch(url,{...options,headers:{'Content-Type':'application/json',...options.headers}});let body;try{body=await response.json();}catch(_){throw new Error(`Unexpected server response (${response.status})`);}if(!response.ok){const error=new Error(typeof body.error==='string'?body.error:body.error?.message||body.message||`HTTP ${response.status}`);error.detail=typeof body.error==='object'?body.error:{message:error.message};throw error;}return body;}
  async function startSimulation() {
    if(state.job){await cancelSimulation();return;}
    if(!project.components.length){failJob({message:'Add components to the canvas first.',stage:'Circuit validation'});return;}
    if(!project.components.some(c=>c.type==='source')){failJob({message:'Add a pressure source to set the circuit pressure reference.',stage:'Circuit validation'});return;}
    document.activeElement?.blur();
    const invalid=$('inspector-content').querySelector(':invalid');if(invalid){invalid.reportValidity();return;}
    clearError();renderAll();state.job='submitting';state.jobRevision=state.revision;setRunning(true);showProgress(0,'Validating circuit…');$('results-panel').classList.remove('collapsed');renderResults();
    try{const data=await api('/api/simulate',{method:'POST',body:JSON.stringify(project)});if(data.result){finishJob(data.result);return;}state.job=data.job_id||data.id;if(!state.job)throw new Error('The server did not return a job ID.');pollJob();}catch(error){failJob(error.detail||{message:error.message});}
  }
  function setRunning(running){$('run-label').textContent=running?'Stop calculation':'Run simulation';$('run-btn').querySelector('.run-icon').textContent=running?'■':'▶';$('run-btn').classList.toggle('running',running);$('engine-select').disabled=running;updateToolbar();if(running){$('result-status').textContent='Calculating';$('result-status').className='status-pill running';}}
  function showProgress(value,message){const progress=Math.max(0,Math.min(1,Number(value)||0));$('progress-value').textContent=`${Math.round(progress*100)}%`;$('progress-bar').style.width=`${progress*100}%`;$('progress-message').textContent=message||'Solving…';$('progress-detail').textContent=project.settings.engine==='comsol'?'COMSOL pipe calculations may take several minutes. You can inspect the schematic or stop the calculation.':'Iterating component flows and node states.';}
  async function pollJob(){if(!state.job||state.job==='submitting')return;const current=state.job;try{const data=await api(`/api/jobs/${encodeURIComponent(current)}`);if(state.job!==current)return;showProgress(data.progress,data.message);if(['complete','completed','done','success'].includes(data.status)){finishJob(data.result);return;}if(data.status==='error'||data.status==='failed'){failJob(typeof data.error==='object'&&data.error?data.error:{message:data.error||data.message||'The solution did not converge.'});return;}if(data.status==='cancelled'){state.job=null;setRunning(false);renderResults();toast('Calculation stopped.');return;}state.pollTimer=setTimeout(pollJob,850);}catch(error){if(state.job===current)failJob(error.detail||{message:error.message});}}
  async function cancelSimulation(){if(!state.job||state.job==='submitting')return;try{$('run-btn').disabled=true;await api(`/api/jobs/${encodeURIComponent(state.job)}/cancel`,{method:'POST',body:'{}'});$('progress-message').textContent='Stopping calculation…';toast('Cancellation requested.');}catch(error){toast(error.message,true);}finally{$('run-btn').disabled=false;}}
  function finishJob(result){state.job=null;setRunning(false);if(!result){failJob('The solver did not return results.');return;}state.result=result;state.resultRevision=state.jobRevision;renderResults();renderScene();renderInspector();toast(state.resultRevision===state.revision?'Simulation complete. Select a component to inspect its profiles.':'Simulation complete. The schematic has changed; results reflect the earlier parameters.');}
  function focusErrorComponent(id) {
    const c=project.components.find(v=>v.id===id);if(!c)return;
    select(id,'params');const v=state.viewport;
    if(c.x-100<v.x||c.x+150>v.x+v.w||c.y-100<v.y||c.y+100>v.y+v.h){v.x=c.x-v.w/2;v.y=c.y-v.h/2;updateViewport();}
    $('inspector-content').scrollTop=0;
  }
  function renderError() {
    const e=state.error,panel=$('simulation-error');panel.hidden=!e;
    $('results-panel').classList.toggle('has-error',!!e);
    if(!e){panel.innerHTML='';return;}
    const components=e.components||[],context=e.context||{};
    const targets=components.map(c=>`<button class="error-component" data-error-component="${esc(c.id)}" ${project.components.some(v=>v.id===c.id)?'':'disabled'} title="Locate ${esc(c.label)} on the schematic">${esc(c.label)}${components.length===1&&context.port?` · Port ${esc(context.port.toUpperCase())}`:''}<span aria-hidden="true"> ↗</span></button>`).join('');
    const details=[e.details||e.message,e.code?`Error code: ${e.code}`:'',...components.map(c=>`Component: ${c.label} (${c.id})`),...Object.entries(context).filter(([k])=>['pressure_a','pressure_b','temperature_a','temperature_b'].includes(k)).map(([k,v])=>`${k}: ${v} ${k.startsWith('pressure')?'Pa':'K'}`)].filter(Boolean).join('\n');
    panel.innerHTML=`<div class="error-heading"><span class="error-mark" aria-hidden="true">!</span><b>${esc(e.stage||'Simulation error')}</b>${targets||'<span class="error-scope">Circuit / solver</span>'}</div><p class="error-message">${esc(e.message)}</p>${e.hint?`<p class="error-hint">${esc(e.hint)}</p>`:''}${e.stale?'<p class="error-hint">This error refers to the submitted schematic. The current circuit has changed.</p>':''}${state.result?'<p class="error-hint">Results shown below are from the previous successful run.</p>':''}<details class="error-details"><summary>Technical details</summary><pre>${esc(details)}</pre></details>`;
  }
  function failJob(detail) {
    clearTimeout(state.pollTimer);const wasJob=!!state.job;state.job=null;
    const e=typeof detail==='string'?{message:detail}:detail;
    e.component_ids=Array.isArray(e.component_ids)?e.component_ids:[];
    e.components=e.components||e.component_ids.map(id=>{const c=project.components.find(v=>v.id===id);return {id,label:c?.label||id};});
    e.stale=wasJob&&state.jobRevision!==state.revision;state.error=e;
    state.validationErrors=e.component_ids;state.errorMessage=e.message;
    setRunning(false);$('results-panel').classList.remove('collapsed');renderResults();renderAll();
    if(e.components.length===1)focusErrorComponent(e.components[0].id);
    document.querySelector('.result-body').scrollTop=0;
    $('simulation-error').focus({preventScroll:true});
  }

  async function confirmAction(title,message){const dialog=$('confirm-dialog');$('confirm-title').textContent=title;$('confirm-message').textContent=message;return new Promise(resolve=>{dialog.addEventListener('close',()=>resolve(dialog.returnValue==='confirm'),{once:true});dialog.showModal();});}
  async function loadExample(){if(state.job){toast('Stop the current calculation first, before replacing the schematic.');return;}if(project.components.length && !await confirmAction('Load example circuit?','The example will replace the current schematic.Export a copy first if needed, Replacement can be undone.'))return;try{const data=await api('/api/example');loadProject(data.project||data,true);toast('Loaded continuous pump Example.');}catch(error){toast(`Failed to load the example: ${error.message}`,true);}}
  function validateProject(data){if(!data||typeof data!=='object'||!Array.isArray(data.components)||!Array.isArray(data.wires))throw new Error('The file must contain components and wires arrays.');if(data.components.length>128||data.wires.length>512)throw new Error('Schematic limit: 128 components, 512 wires.');const ids=new Set();for(const c of data.components){if(!c||!Object.hasOwn(types,c.type))throw new Error(`Unsupported component type: ${c?.type}`);if(!c.id||typeof c.id!=='string'||ids.has(c.id))throw new Error('Component ID missing or duplicated.');ids.add(c.id);if(!Number.isFinite(c.x)||!Number.isFinite(c.y))throw new Error('Invalid component position.');if(!c.params||typeof c.params!=='object'||Array.isArray(c.params))throw new Error('Invalid component parameters.');c.label=String(c.label||types[c.type].prefix).slice(0,80);c.rotation=Number.isFinite(c.rotation)?Math.round(c.rotation/90)*90:0;if(c.flipX!==undefined&&typeof c.flipX!=='boolean')throw new Error('Component flipX must be true or false.');c.flipX=c.flipX??false;if(c.type==='jt'){c.params.law??='orifice';c.params.impedance_mode??='geometry';c.params.material??='';}c.params={...clone(types[c.type].params),...c.params};
      for(const [key,value] of Object.entries(types[c.type].params)){if(typeof value==='number'&&!Number.isFinite(c.params[key]))throw new Error(`${c.label} of ${key} must be a number.`);}
      if(c.type==='pump')for(const key of ['forward_speed_curve','compression_curve']){if(!validPumpCurve(c.params[key],key))throw new Error(`${c.label}: invalid ${key}. Use positive pressure/value pairs with strictly increasing pressures; compression must be at least one.`);}
      if(c.type==='pipe' && c.params.wall_boundary==='prescribed'){
        const rows=c.params.wall_profile;if(!Array.isArray(rows)||!rows.length||rows.length>32)throw new Error('A prescribed wall profile must contain 1–32 temperature segments.');
        for(const row of rows){if(!row||!['length_m','T_start_K','T_end_K'].every(key=>Number.isFinite(row[key])&&row[key]>0))throw new Error('Wall-profile segment lengths and temperatures must be positive.');}
      }
    }
    const wireIds=new Set();for(const w of data.wires){if(!w.id)w.id=uid('w');if(wireIds.has(w.id)||ids.has(w.id))throw new Error('Wire ID Duplicate.');wireIds.add(w.id);for(const ep of [w.from,w.to]){if(!ep||!ids.has(ep.component)||!['a','b'].includes(ep.port))throw new Error('A wire references an invalid port.');const c=data.components.find(c=>c.id===ep.component);if(c.type==='junction'&&ep.port!=='a')throw new Error('A junction has only port A.');}}
    data.settings={...data.settings,engine:['comsol','reduced','two_phase'].includes(data.settings?.engine)?data.settings.engine:'comsol'};data.name=String(data.name||'Import circuit').slice(0,80);data.version=1;return data;
  }
  function loadProject(data,history=true){const checked=validateProject(clone(data));if(history)commit();project=checked;state.selected=null;state.pending=null;cancelConnection();resetResult();changed();syncControls();requestAnimationFrame(fit);}
  function download(filename,data){const blob=new Blob([JSON.stringify(data,null,2)],{type:'application/json'});const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}

  function bindCanvas(){
    svg.addEventListener('pointerdown',event=>{
      state.canvasFocused=true;svg.focus();const p=worldPoint(event);
      if(event.button===1 || (event.button===0&&(state.tool==='pan'||state.space))){event.preventDefault();state.pointer={kind:'pan',clientX:event.clientX,clientY:event.clientY,view:clone(state.viewport)};svg.setPointerCapture(event.pointerId);$('canvas-wrap').classList.add('panning');return;}
      if(event.button!==0)return;
      const port=event.target.closest('[data-port]');if(port){event.stopPropagation();connectPort(port.dataset.component,port.dataset.port);return;}
      const component=event.target.closest('[data-id]');if(component){const c=project.components.find(v=>v.id===component.dataset.id);select(c.id);state.pointer={kind:'drag',id:c.id,startX:p.x,startY:p.y,x:c.x,y:c.y,moved:false,snapshot:JSON.stringify(project)};svg.setPointerCapture(event.pointerId);return;}
      const wire=event.target.closest('[data-wire]');if(wire){select(wire.dataset.wire);return;}
      if(state.pending){cancelConnection();return;}select(null);
    });
    svg.addEventListener('pointermove',event=>{
      const p=worldPoint(event);
      if(state.pending){const a=endpoint(state.pending);if(a){$('pending-wire').removeAttribute('hidden');$('pending-wire').setAttribute('d',wirePath(a,p,state.pending));}}
      const move=state.pointer;if(!move)return;
      if(move.kind==='pan'){const rect=svg.getBoundingClientRect();const scale=move.view.w/rect.width;state.viewport.x=move.view.x-(event.clientX-move.clientX)*scale;state.viewport.y=move.view.y-(event.clientY-move.clientY)*scale;updateViewport();}
      else {const c=project.components.find(v=>v.id===move.id);if(!c)return;const dx=p.x-move.startX,dy=p.y-move.startY;if(Math.abs(dx)+Math.abs(dy)>4){move.moved=true;c.x=snap(move.x+dx);c.y=snap(move.y+dy);renderScene();}}
    });
    function endPointer(event){const move=state.pointer;if(!move)return;state.pointer=null;$('canvas-wrap').classList.remove('panning');if(svg.hasPointerCapture(event.pointerId))svg.releasePointerCapture(event.pointerId);if(move.kind==='drag'&&move.moved){state.undo.push(move.snapshot);if(state.undo.length>60)state.undo.shift();state.redo=[];changed({render:false});}}
    svg.addEventListener('pointerup',endPointer);svg.addEventListener('pointercancel',endPointer);
    svg.addEventListener('wheel',event=>{event.preventDefault();zoom(event.deltaY>0?1.1:1/1.1,worldPoint(event));},{passive:false});
    svg.addEventListener('contextmenu',event=>event.preventDefault());
    $('canvas-wrap').addEventListener('dragover',event=>{if(event.dataTransfer.types.includes('application/x-helium-component')){event.preventDefault();event.dataTransfer.dropEffect='copy';}});
    $('canvas-wrap').addEventListener('drop',event=>{const type=event.dataTransfer.getData('application/x-helium-component');if(types[type]){event.preventDefault();addComponent(type,worldPoint(event));}});
    document.addEventListener('keydown',event=>{
      const tag=document.activeElement?.tagName;if(['INPUT','TEXTAREA','SELECT'].includes(tag)||document.activeElement?.isContentEditable||$('confirm-dialog').open)return;
      if(!state.canvasFocused)return;const key=event.key.toLowerCase();
      if(event.code==='Space'){event.preventDefault();state.space=true;$('canvas-wrap').classList.add('pan');}
      else if(key==='escape'){cancelConnection();select(null);}
      else if((event.ctrlKey||event.metaKey)&&key==='z'){event.preventDefault();event.shiftKey?redo():undo();}
      else if((event.ctrlKey||event.metaKey)&&key==='y'){event.preventDefault();redo();}
      else if(key==='delete'||key==='backspace'){event.preventDefault();removeSelected();}
      else if(!event.ctrlKey&&!event.metaKey&&!event.altKey){if(key==='r')rotateSelected();else if(key==='x'){event.preventDefault();flipSelected();}else if(key==='v')setTool('select');else if(key==='h')setTool('pan');else if(key==='f')fit();}
    });
    document.addEventListener('keyup',event=>{if(event.code==='Space'){state.space=false;$('canvas-wrap').classList.toggle('pan',state.tool==='pan');}});
    window.addEventListener('blur',()=>{state.space=false;$('canvas-wrap').classList.toggle('pan',state.tool==='pan');});
    // Selection replaces SVG nodes during this event; use its original path.
    document.addEventListener('pointerdown',event=>{if(!event.composedPath().some(el=>el instanceof Element&&el.matches('.canvas-wrap,.canvas-toolbar,.palette')))state.canvasFocused=false;});
    let resizeTimer;new ResizeObserver(()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(()=>{const rect=svg.getBoundingClientRect();if(rect.height>0){state.viewport.h=state.viewport.w*rect.height/rect.width;updateViewport();}},70);}).observe($('canvas-wrap'));
  }
  function bindControls(){
    $('palette').innerHTML=Object.entries(types).map(([type,t])=>`<button class="palette-item" draggable="true" data-add="${type}" data-testid="add-${type}" title="Add${t.name}" aria-label="Add${t.name}"><svg viewBox="0 0 40 28">${icons[type]}</svg><span><b>${t.name}</b><small>${t.en}</small></span></button>`).join('');
    $('palette').querySelectorAll('[data-add]').forEach(button=>{button.onclick=()=>addComponent(button.dataset.add);button.ondragstart=event=>{event.dataTransfer.setData('application/x-helium-component',button.dataset.add);event.dataTransfer.effectAllowed='copy';};});
    $('select-tool').onclick=()=>setTool('select');$('pan-tool').onclick=()=>setTool('pan');$('undo-btn').onclick=undo;$('redo-btn').onclick=redo;$('rotate-btn').onclick=rotateSelected;$('flip-btn').onclick=flipSelected;$('delete-btn').onclick=removeSelected;
    $('zoom-in').onclick=()=>zoom(1/1.2);$('zoom-out').onclick=()=>zoom(1.2);$('fit-btn').onclick=fit;$('zoom-reset').onclick=fit;
    $('project-name').onchange=event=>{const name=event.target.value.trim();if(name&&name!==project.name){commit();project.name=name;changed({render:false});}else event.target.value=project.name;};
    $('engine-select').onchange=event=>{commit();project.settings.engine=event.target.value;changed();};
    $('params-tab').onclick=()=>{state.tab='params';renderInspector();};$('outputs-tab').onclick=()=>{state.tab='outputs';renderInspector();};
    $('run-btn').onclick=startSimulation;$('example-btn').onclick=loadExample;$('empty-example-btn').onclick=loadExample;
    $('clear-btn').onclick=async()=>{if(state.job){toast('Stop the current calculation first, before replacing the schematic.');return;}if(project.components.length&&!await confirmAction('New schematic?','The current schematic will be cleared.Export a copy first if needed, Clearing can be undone.'))return;commit();project=freshProject();state.selected=null;cancelConnection();resetResult();changed();syncControls();fit();};
    $('export-btn').onclick=()=>download(`${project.name.replace(/[<>:"/\\|?*]/g,'_')||'fluid_circuit'}.json`,project);
    $('import-btn').onclick=()=>{if(state.job){toast('Stop the current calculation first, before importing a schematic.');return;}$('import-file').click();};
    $('import-file').onchange=async event=>{const file=event.target.files[0];if(!file)return;try{if(file.size>1024*1024)throw new Error('File exceeds 1 MB limits.');const data=JSON.parse(await file.text());validateProject(clone(data));if(project.components.length&&!await confirmAction('Import circuit file?','Importing will replace the current schematic, Replacement can be undone.'))return;loadProject(data);toast('Circuit imported.');}catch(error){toast(`Import failed: ${error.message}`,true);}finally{event.target.value='';}};
    $('simulation-error').onclick=event=>{const button=event.target.closest('[data-error-component]');if(button)focusErrorComponent(button.dataset.errorComponent);};
    $('results-toggle').onclick=()=>{$('results-panel').classList.toggle('collapsed');};
    $('color-by').onchange=event=>{state.color=event.target.value;renderScene();};
    $('download-results').onclick=()=>{if(state.result)download(`${project.name}_results.json`,{project:clone(project),result:state.result,stale:state.resultRevision!==state.revision});};
    $('result-table').onclick=event=>{const row=event.target.closest('[data-result-id]');if(row)select(row.dataset.resultId,'outputs');};
  }
  async function initialize(){bindControls();bindCanvas();try{const saved=localStorage.getItem(STORE);if(saved)project=validateProject(JSON.parse(saved));}catch(_){toast('The local draft could not be read, New schematic created.',true);}syncControls();renderAll();renderResults();requestAnimationFrame(fit);
    try{const status=await api('/api/status');state.status=status;$('service-dot').classList.add('online');$('service-label').textContent='Local solver connected';}catch(_){try{await api('/api/example');$('service-dot').classList.add('online');$('service-label').textContent='Local solver connected';}catch(error){$('service-label').textContent='Solver disconnected';}}
  }
  // A read-only inspection surface makes browser verification independent of internal UI state.
  window.HeliumStudio={getProject:()=>clone(project),getResult:()=>clone(state.result),getSelected:()=>state.selected};
  initialize();
})();
