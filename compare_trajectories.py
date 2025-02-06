import numpy as np

import argparse
import matplotlib.pyplot as plt
from particle_in_mag import track_particle
from geant4 import run as run_g4


def plot_trajectories_rk(ax, initial_conditions, num_steps=8000, mag_field='toy'):
    trajectories = []
    for conditions in initial_conditions:
        x, y, z, px, py, pz, charge = conditions
        trajectory = track_particle(x, y, z, px, py, pz, charge, num_steps, mag_field=mag_field)
        trajectories.append(trajectory)

    colors = ['lightblue', 'lightcoral', 'lightgreen', 'gray']
    labels = ['Q= -1, PYTHON', 'Q=1, PYTHON', 'Muon3 PYTHON', 'Muon4 PYTHON']
    for i, trajectory in enumerate(trajectories):
        positions = np.array(trajectory)[:, :3]
        ax.plot(positions[:, 2], positions[:, 0], positions[:, 1], color=colors[i], label=labels[i], linestyle = 'dashed')
        final_position = trajectory[-1][:3]
        final_momentum = trajectory[-1][3:]
        print(f"{labels[i]} final position: {final_position}, final momentum: {final_momentum}")

    return trajectories

def plot_trajectories_g4(ax, initial_conditions, mag_field='toy'):
    muon_data = run_g4(initial_conditions, mag_field)
    colors = ['b', 'r', 'g', 'k']
    labels = ['Q= -1, GEANT4', 'Q=1, GEANT4', 'Muon3 GEANT4', 'Muon4 GEANT4']
    for i, muon in enumerate(muon_data):
        positions = np.array([muon['x'], muon['y'], muon['z']]).T
        ax.plot(positions[:, 2], positions[:, 0], positions[:, 1], color=colors[i % len(colors)], label=labels[i])
        final_position = positions[-1]
        final_momentum = np.array([muon['px'][-1], muon['py'][-1], muon['pz'][-1]])
        print(f"{labels[i]} final position: {final_position}, final momentum: {final_momentum}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag_type", type=str, default='toy')
    args = parser.parse_args()

    data = np.array([
        [0.0, 0.0, -1, 0.0, 0.0, 20.0, -1.0], 
        [0.0, 0.0, -1, 0.0, 0.0, 20.0, 1.0],    
        [-0.1, 0.1, -1, 1.0, 1.0, 20.0, -1.0], 
        [-2., -1., 10., 4, 0.7, 15.0, 1.0],   
    ])

    if args.mag_type == 'uniform_map':
        mag_type = '/home/hep/lprate/projects/MuonsAndMatter/data/outputs/uniform_fields.pkl'
    elif args.mag_type == 'map':
        mag_type = '/home/hep/lprate/projects/MuonsAndMatter/data/outputs/fields.pkl'
    else:
        mag_type = args.mag_type
    


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    trajectories_rk =  plot_trajectories_rk(ax,data, 8000, mag_type)
    trajectories_g4 = plot_trajectories_g4(ax, data, mag_type)

    ax.set_xlabel('Z (m)')
    ax.set_ylabel('X (m)')
    ax.set_zlabel('Y (m)')
    ax.set_title('3D Trajectories of Muons')

    ax.set_xlim(100, -5)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    ax.view_init(elev=60, azim=110)

    ax.legend(loc='upper left', bbox_to_anchor=(-0.3, 1))
    plt.savefig('trajectory_geant.png')
    plt.show()
    
    