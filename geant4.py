import json
import numpy as np
from muon_slabs import initialize, simulate_muon, collect
from time import time
import pickle
from mag_fields import UniformMagneticField

def get_design(d:dict):
    detector = {
         "worldPositionX": 0, "worldPositionY": 0, "worldPositionZ": 0, "worldSizeX": 21, "worldSizeY": 21,
         "worldSizeZ": 200,
        "type": 4,
        "limits" : {
            "max_step_length": 0.05,
            "minimum_kinetic_energy": 0.1,
        },
        "global_field_map": d,
    }
    return detector

def initialize_geant4(detector, seed=None):
    B = np.array(detector['global_field_map'].pop('B'))
    if seed is None:
        seeds = (np.random.randint(256), np.random.randint(256), np.random.randint(256), np.random.randint(256))
    else:
        seeds = (seed, seed, seed, seed)
    output_data = initialize(*seeds, json.dumps(detector), B)
    return output_data

def simulate_muons(muons):
    muon_data = []
    for muon in muons:
        x, y, z,px,py,pz, charge = muon[:7]
        simulate_muon(px, py, pz, int(charge), x, y, z)
        data = collect()
        muon_data.append(data)
    return muon_data

def get_field_dict(file_name=None):
    with open(file_name, 'rb') as f:
        fields = pickle.load(f)
    points,B = fields['points'], fields['B'].astype(np.float32)
    d_space = (points[:, 0].max().item(), points[:, 1].max().item(), (points[:, 2].min().item(), points[:, 2].max().item()))
    resol = (np.diff(np.unique(points[:, 0]))[0].item(), np.diff(np.unique(points[:, 1]))[0].item(), np.diff(np.unique(points[:, 2]))[0].item())
    return {'B': B,
            'range_x': [0.,d_space[0], resol[0]],
            'range_y': [0.,d_space[1], resol[1]],
            'range_z': [d_space[2][0],d_space[2][1], resol[2]]}

def run(data,mag_type:str):
    if mag_type == 'toy':
        field_map = {'B': []}
    elif mag_type == 'uniform':
        field_map = {'B': UniformMagneticField.get_magnetic_field(0, 0, 0)}
    else:
        field_map = get_field_dict(file_name=mag_type)
        
    detector = get_design(field_map)
    detector["store_primary"] = True
    detector["store_all"] = False
    
    t1_init = time()
    output_data = initialize_geant4(detector, 10)
    print('Time to initialize', time() - t1_init)
    muon_data = simulate_muons(data)
    print(f'Simulated {len(data)} muons')
    return muon_data