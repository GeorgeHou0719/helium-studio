"""Build the user's stationary He4/304 pipe in COMSOL 6.1 through MPh.

Run from workspace root: python -m models.comsol_pipe.build_model
The resulting continuum model is a formal baseline at the specified vacuum
outlet: Kn is about 8, outside continuum applicability.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import time

import numpy as np
from jpype import JInt
from scipy.integrate import cumulative_trapezoid

from models.one_k_pot.he4_viscosity import gas_viscosity
from .steel_properties import comsol_conductivity_expression, thermal_conductivity_w_m_k

HERE = Path(__file__).resolve().parent
from models.config import COMSOL_EXECUTABLE as COMSOL, COMSOL_TEMPLATE


def interpolation(j, tag, rows, argunit, fununit, label):
    f = j.func().create(tag, 'Interpolation')
    f.label(label)
    f.set('funcname', tag)
    f.set('source', 'table')
    f.set('table', [[f'{x:.16g}', f'{y:.16g}'] for x, y in rows])
    f.set('argunit', argunit)
    f.set('fununit', fununit)
    f.set('interp', 'piecewisecubic')
    f.set('extrap', 'const')
    return f


def box(co, tag, dim, bounds, label):
    s = co.selection().create(tag, 'Box')
    s.label(label)
    s.set('entitydim', JInt(dim))
    s.set('condition', 'inside')
    for prop, value in zip(('xmin', 'xmax', 'ymin', 'ymax'), bounds):
        s.set(prop, value)
    return s


def build(client, nz=640, nr=32, conductivity='nist'):
    from .helium_properties import gas_thermal_conductivity
    m = client.create('Steady helium flow in insulated SS304 pipe')
    j = m.java
    j.modelPath(str(HERE))
    j.label('helium_ss304_pipe.mph')
    j.comments('User-specified stationary continuum laminar baseline. L=288 mm, Do=6.35 mm, wall=0.762 mm; absolute pressures 15 and 0.5 Pa; end temperatures 1 and 293.15 K. Outer wall insulated. Outlet Kn~8: continuum/no-slip predictions are not validated physical transport. See companion README.')
    params = {
        'L': ('288[mm]', 'Pipe length'),
        'Do': ('6.35[mm]', 'Outer diameter'),
        'wall_w': ('0.762[mm]', 'Wall thickness w (w is reserved for axial velocity)'),
        'Ri': ('Do/2-wall_w', 'Inner radius'),
        'Di': ('2*Ri', 'Bore diameter'),
        'P1': ('15[Pa]', 'Cold inlet ABSOLUTE pressure'),
        'P2': ('0.5[Pa]', 'Warm outlet ABSOLUTE pressure'),
        'T1': ('1[K]', 'Cold end gas and wall temperature'),
        'T2': ('293.15[K]', 'Warm end gas and wall temperature'),
        'MHe': ('0.004002602[kg/mol]', 'He4 molar mass, inherited Python model'),
        'RHe': ('8.314462618[J/(mol*K)]/MHe', 'Specific ideal-gas constant'),
        'cpHe': ('2.5*RHe', 'Dilute monatomic gas heat capacity'),
        'load': ('1', 'Continuation factor: 1 gives requested inlet pressure'),
        'P1run': ('P2+load*(P1-P2)', 'Continuation inlet pressure'),
        'epssel': ('1e-6[m]', 'Selection tolerance (larger than rendering roundoff, smaller than wall)'),
    }
    for key, (value, desc) in params.items():
        j.param().set(key, value, desc)
    temperatures = np.unique(np.r_[np.geomspace(1, 10, 250), np.linspace(10, 350, 1000)])
    interpolation(j, 'muHe', [(t, gas_viscosity(float(t))) for t in temperatures], 'K', 'Pa*s', 'He4 viscosity: existing one_k_pot NIST Hermite model')
    kt = temperatures if conductivity == 'nist' else np.unique(np.r_[temperatures[temperatures<300],300.])
    interpolation(j, 'kHe', [(t, gas_thermal_conductivity(float(t), source=conductivity)) for t in kt], 'K', 'W/(m*K)', f'He4 thermal conductivity: {conductivity}')
    f = j.func().create('kSS', 'Analytic')
    f.label('304 stainless: NIST conductivity fit, equation range 1-300 K')
    f.set('funcname', 'kSS')
    f.set('args', 'theta')
    # Positive clamp is a Newton-iteration guard, inactive for the solved range.
    f.set('expr', comsol_conductivity_expression('max(theta,1)').removesuffix('*1[W/(m*K)]'))
    f.set('argunit', 'K')
    f.set('fununit', 'W/(m*K)')

    # Conductive initial guess only; the final temperature is solved by COMSOL.
    tt = np.geomspace(1, 293.15, 8000)
    integ = cumulative_trapezoid(thermal_conductivity_w_m_k(tt), tt, initial=0)
    xx = integ/integ[-1]
    interpolation(j, 'Tguess', zip(xx, tt), '1', 'K', 'Initial guess from integrated steel conductivity')
    resistance = cumulative_trapezoid(np.array([gas_viscosity(float(t)) for t in tt])*tt, xx, initial=0)
    pg = np.sqrt(15**2-(15**2-.5**2)*resistance/resistance[-1])
    interpolation(j, 'pguess', zip(xx, pg), '1', 'Pa', 'Initial guess only: pressure-squared Poiseuille integral')
    mdot = np.pi*.004826**4*(225-.25)/(256*(8.314462618/.004002602)*.288*resistance[-1])
    j.param().set('mdot_guess', f'{mdot:.16g}[kg/s]', 'Initialization only, final mass flow is solved')
    j.component().create('comp1', True)
    co = j.component('comp1')
    g = co.geom().create('geom1', 2)
    g.axisymmetric(True)
    g.lengthUnit('m')
    for tag, pos, size, label in (
        ('gas', ['0','0'], ['Ri','L'], 'Helium gas bore'),
        ('steel', ['Ri','0'], ['wall_w','L'], '304 stainless wall')):
        rect = g.create(tag, 'Rectangle')
        rect.label(label)
        rect.set('pos', pos)
        rect.set('size', size)
        rect.set('selresult', True)
    g.run()
    selections = [
        ('inlet', 1, ['-epssel','Ri+epssel','-epssel','epssel'], 'Gas inlet, z=0'),
        ('outlet', 1, ['-epssel','Ri+epssel','L-epssel','L+epssel'], 'Gas outlet, z=L'),
        ('cold', 1, ['-epssel','Do/2+epssel','-epssel','epssel'], 'Cold ends: gas and steel'),
        ('warm', 1, ['-epssel','Do/2+epssel','L-epssel','L+epssel'], 'Warm ends: gas and steel'),
        ('outer', 1, ['Do/2-epssel','Do/2+epssel','-epssel','L+epssel'], 'Insulated outer wall'),
        ('interface', 1, ['Ri-epssel','Ri+epssel','-epssel','L+epssel'], 'Gas/steel shared interface'),
        ('axis', 1, ['-epssel','epssel','-epssel','L+epssel'], 'Symmetry axis'),
        ('steelends', 1, ['Ri-epssel','Do/2+epssel','-epssel','epssel'], 'Steel cold end'),
        ('steelwarm', 1, ['Ri-epssel','Do/2+epssel','L-epssel','L+epssel'], 'Steel warm end'),
    ]
    expected = {'inlet':[2], 'outlet':[3], 'cold':[2,5], 'warm':[3,6], 'outer':[7], 'interface':[4], 'axis':[1], 'steelends':[5], 'steelwarm':[6]}
    for tag, dim, bounds, label in selections:
        s=box(co, tag, dim, bounds, label)
        print('SELECTION', tag, list(s.entities(dim)), flush=True)
        assert list(s.entities(dim)) == expected[tag], f'Unexpected selection: {tag}'

    mat = co.material().create('he', 'Common')
    mat.label('Dilute helium-4 gas: inherited viscosity and ideal gas; NIST conductivity')
    mat.selection().named('geom1_gas_dom')
    mat.materialType('nonSolid')
    for prop, val in {'density':'p/(RHe*T)', 'dynamicviscosity':'muHe(T)', 'thermalconductivity':'kHe(T)', 'heatcapacity':'cpHe', 'ratioofspecificheat':'5/3', 'molarmass':'MHe', 'thermalexpansioncoefficient':'1/T', 'soundspeed':'sqrt(5/3*RHe*T)'}.items():
        mat.propertyGroup('def').set(prop, val)
    mat = co.material().create('ss', 'Common')
    mat.label('SS304: NIST k(T); room-temperature rho/Cp unused in stationary solid')
    mat.selection().named('geom1_steel_dom')
    for prop, val in {'density':'7900[kg/m^3]', 'heatcapacity':'500[J/(kg*K)]', 'thermalconductivity':'kSS(T)'}.items():
        mat.propertyGroup('def').set(prop, val)

    spf=co.physics().create('spf','LaminarFlow','geom1')
    spf.label('Laminar helium flow - absolute pressure')
    spf.selection().named('geom1_gas_dom')
    spf.prop('PhysicalModelProperty').set('Compressibility','CompressibleMALT03')
    spf.prop('PhysicalModelProperty').set('pref','0[Pa]')
    spf.prop('ShapeProperty').set('order_fluid','2')
    spf.feature('fp1').set('rho_mat','userdef')
    spf.feature('fp1').set('rho','p/(RHe*T)')
    spf.feature('fp1').set('mu_mat','userdef')
    spf.feature('fp1').set('mu','muHe(T)')
    spf.feature('init1').set('p_init','pguess(z/L)')
    spf.feature('init1').set('u_init',['0','0','0'])
    inlet=spf.create('in1','InletBoundary',1)
    inlet.label('Cold inlet - prescribed mean pressure P1')
    inlet.selection().named('inlet')
    inlet.set('BoundaryCondition','FullyDevelopedFlow')
    inlet.set('FullyDevelopedFlowOption','p0av')
    inlet.set('p0avfdf','P1run')
    inlet.set('SuppressBackflow',False)
    out=spf.create('out1','OutletBoundary',1)
    out.label('Warm outlet - prescribed mean pressure P2')
    out.selection().named('outlet')
    out.set('BoundaryCondition','FullyDevelopedFlow')
    out.set('FullyDevelopedFlowOption','p0av')
    out.set('p0avfdf','P2')
    out.set('SuppressBackflow',False)
    ht=co.physics().create('ht','HeatTransferInSolidsAndFluids','geom1')
    ht.label('Conjugate heat transfer - helium and SS304')
    print('Flow physics configured',flush=True)
    # Default solid node is noneditable; the selected Fluid node overrides it.
    ht.feature('fluid1').selection().named('geom1_gas_dom')
    ht.feature('fluid1').set('fluidType','idealGas')
    ht.feature('fluid1').set('Cp_mat','userdef')
    ht.feature('fluid1').set('Cp','cpHe')
    ht.feature('fluid1').set('Rs_mat','userdef')
    ht.feature('fluid1').set('Rs','RHe')
    ht.feature('init1').set('Tinit','Tguess(z/L)')
    for tag, selection, temp in [('temp1','cold','T1'),('temp2','warm','T2')]:
        f=ht.create(tag,'TemperatureBoundary',1)
        f.label(f'{selection.title()} end temperature - gas and wall')
        f.selection().named(selection)
        f.set('T0',temp)
    n=co.multiphysics().create('nitf1','NonIsothermalFlow',2)
    # Coupling domain is automatically the intersection of fluid selections.
    n.set('Fluid_physics','spf'); n.set('Heat_physics','ht')
    n.set('SpecifyDensity','FromFluidFlowInterface')
    n.set('includeViscousDissipation', True)
    print('Heat physics and coupling configured',flush=True)
    v=co.variable().create('diagnostics')
    v.selection().named('geom1_gas_dom')
    for key, expr in {'rhoHe':'p/(RHe*T)', 'lambdaHe':'muHe(T)/p*sqrt(pi*RHe*T/2)', 'Kn':'lambdaHe/Di', 'Mach':'spf.U/sqrt(5/3*RHe*T)', 'Re_local':'rhoHe*spf.U*Di/muHe(T)', 'massflux':'rhoHe*w'}.items():
        v.set(key,expr)
    for tag, sel in [('intin','inlet'),('intout','outlet'),('intcold','cold'),('intwarm','warm'),('intsc','steelends'),('intsw','steelwarm')]:
        op=co.cpl().create(tag,'Integration')
        op.selection().geom('geom1',1)
        op.selection().named(sel)
    # Custom coupling integration uses the meridional measure. Expressions in
    # postprocessing include 2*pi*r explicitly and verify the resulting area.
    mesh=co.mesh().create('mesh1')
    mp=mesh.create('map1','Map')
    ax=mp.create('axial','Distribution')
    ax.selection().set([1,4,7])
    ax.set('type','predefined'); ax.set('elemcount',JInt(nz))
    ax.set('elemratio',10000.0); ax.set('method','geometric'); ax.set('symmetric',True)
    for tag, edges, count in [('radialgas',[2,3],nr),('radialsteel',[5,6],max(4,nr//4))]:
        d=mp.create(tag,'Distribution'); d.selection().set(edges); d.set('numelem',JInt(count))
    mesh.run()
    print('MESH complete',nz,nr, flush=True)
    st=j.study().create('std1'); st.label('Stationary conjugate laminar baseline')
    stat=st.create('stat','Stationary')
    j.sol().create('sol1')
    sol=j.sol('sol1'); sol.study('std1'); sol.attach('std1')
    sol.create('st1','StudyStep'); sol.feature('st1').set('study','std1'); sol.feature('st1').set('studystep','stat')
    sol.create('v1','Variables')
    so=sol.create('s1','Stationary')
    so.set('stol',1e-6)
    so.create('fc1','FullyCoupled'); so.feature('fc1').set('maxiter',JInt(80))
    so.create('d1','Direct'); so.feature('d1').set('linsolver','pardiso')
    so.feature('fc1').set('linsolver','d1')
    so.feature().remove('fcDef')
    for tag, expr, label in [('pgT','T','Temperature (K)'),('pgU','spf.U','Gas velocity (m/s)'),('pgP','p','Absolute pressure (Pa)'),('pgKn','Kn','Knudsen number: continuum validity')]:
        pg=j.result().create(tag,'PlotGroup2D'); pg.label(label)
        surf=pg.create('surf1','Surface'); surf.set('expr',expr)
    return m


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-only',action='store_true')
    parser.add_argument('--nz',type=int,default=640)
    parser.add_argument('--nr',type=int,default=32)
    parser.add_argument('--conductivity',choices=['nist','legacy'],default='nist')
    parser.add_argument('--output',type=Path,default=COMSOL_TEMPLATE)
    args=parser.parse_args()
    if not COMSOL.is_file():
        parser.error('Set COMSOL_EXECUTABLE to the COMSOL 6.1 executable before building the template.')
    os.environ['PATH']=str(COMSOL.parent)+os.pathsep+os.environ['PATH']
    import mph
    start=time.time()
    client=mph.Client(cores=4,version='6.1')
    print('COMSOL ready',flush=True)
    m=build(client,args.nz,args.nr,args.conductivity)
    output=args.output.resolve(); output.parent.mkdir(parents=True,exist_ok=True)
    m.save(str(output))
    print('SAVED configured model',output,flush=True)
    if not args.build_only:
        print('SOLVING',flush=True)
        m.java.sol('sol1').runAll()
        m.save(str(output))
        print('SOLVED',flush=True)
        for expr in ['comp1.intin(2*pi*r*comp1.rhoHe*w)','comp1.intout(2*pi*r*comp1.rhoHe*w)']:
            print(expr,m.evaluate(expr),flush=True)
        from .inspect_results import inspect
        inspect(m, output.parent/'results')
        m.save(str(output))
    print('ELAPSED',time.time()-start,flush=True)


if __name__=='__main__':
    main()
