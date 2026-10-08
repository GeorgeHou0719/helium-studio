"""Editable examples and component defaults; all internal values use SI."""
from copy import deepcopy
from .pumps import FORWARD_BACKFLOW_DEFAULTS

WALL_PROFILE=[
    {'length_m':.3,'T_start_K':293.15,'T_end_K':293.15},
    {'length_m':.3,'T_start_K':293.15,'T_end_K':35.},
    {'length_m':.3,'T_start_K':35.,'T_end_K':35.},
    {'length_m':.3,'T_start_K':35.,'T_end_K':4.},
    {'length_m':.3,'T_start_K':4.,'T_end_K':4.},
]
DEFAULTS={
    'pipe':{'length_m':.288,'outer_diameter_m':.00635,'wall_thickness_m':.000762,
            'T_a_K':293.15,'T_b_K':293.15,'wall_boundary':'adiabatic'},
    'pump':{'mode':'constant_speed','speed_m3_s':.004,'expression':'Pin*S/(R*Tin)',
            'outlet_temperature_K':293.15,**deepcopy(FORWARD_BACKFLOW_DEFAULTS)},
    'source':{'mode':'constant','pressure_Pa':15.57,'temperature_K':1.,'expression':'15.57','heat_load_W':0.},
    'jt':{'law':'capillary','material':'CuNi','diameter_m':1e-4,'length_m':.1,
          'impedance_mode':'measured','impedance_m_inv3':3.2e16,'discharge_coefficient':.7},
    'junction':{'temperature_K':293.15},
}


def continuous_pump_example():
    def component(cid,kind,label,x,y,rotation=0,**params):
        return dict(id=cid,type=kind,label=label,x=x,y=y,rotation=rotation,
                    params={**deepcopy(DEFAULTS[kind]),**params})
    parts=[
        component('pot','source','Saturated vapor pot',180,300,270,mode='saturation'),
        component('pipe1','pipe','Suction pipe',440,160,length_m=.5,outer_diameter_m=.0635,
                  wall_thickness_m=.00762,T_a_K=1.,T_b_K=293.15),
        component('pump1','pump','Circulation pump',740,160,mode='forward_backflow',
                  forward_speed_m3_s=.240,compression_ratio=1e8,outlet_temperature_K=293.15),
        component('pipe2','pipe','Staged precooling pipe',740,440,180,length_m=1.5,T_a_K=293.15,T_b_K=4.,
                  outer_diameter_m=.0635,wall_thickness_m=.00762,
                  wall_boundary='prescribed',wall_profile=deepcopy(WALL_PROFILE)),
        component('jt1','jt','CuNi J–T capillary',440,440,180),
    ]
    pairs=[('pot','b','pipe1','a'),('pipe1','b','pump1','a'),('pump1','b','pipe2','a'),
           ('pipe2','b','jt1','a'),('jt1','b','pot','a')]
    wires=[{'id':f'w{i+1}','from':{'component':a,'port':ap},'to':{'component':b,'port':bp}}
           for i,(a,ap,b,bp) in enumerate(pairs)]
    return {'version':1,'name':'Continuous pump · Staged precooling circuit',
            'settings':{'engine':'comsol'},'components':parts,'wires':wires,
            'assumptions':[
                'User specified: Precooler outer-wall total length 1.5 m; five 300 mm segments: 293.15 K, 293.15→35 K, 35 K, 35→4 K, 4 K. J–T inlet boundary: 4 K.',
                'User specified: Both 304 stainless-steel pipes have outer diameter 63.5 mm and wall thickness 7.62 mm; suction length 500 mm, adiabatic outer wall, endpoints 1 K and 293.15 K.',
                'User specified: Saturated vapor pot at 1 K; forward–backflow pump with nominal forward speed 240 L/s, compression ratio 1e8, and discharge temperature 293.15 K.',
                'User-selected literature restriction: CuNi capillary, bore 0.10 mm, length 100 mm, measured impedance z=3.2e16 m^-3. This is not a measurement of the user apparatus.',
                'No conduction is exchanged between components; gas enthalpy is retained in the pot balance. Two-phase hydraulics require the Two-phase pot (1D) engine.',
            ]}
