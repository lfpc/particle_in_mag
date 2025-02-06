import numpy as np
import pickle

MUON_MASS = 0.1056583755  # GeV/c²
c = 299792458.0  # Speed of light in m/s
e_charge = 1.602176634e-19  # Elementary charge in Coulombs
# Conversion factor: 1 GeV/c in SI momentum units (kg·m/s)
GeV_over_c_to_SI = 1.602176634e-10 / c  # ≈ 5.3443e-19 kg·m/s


from mag_fields import ToyMagneticField, UniformMagneticField, CustomMagneticField
    
def rk4_step(x: float, y: float, z: float, 
             px: float, py: float, pz: float, 
             charge: float, step_length: float,
             mag_field):
    """
    Performs one RK4 integration step with fixed spatial step size.
    
    Args:
        x, y, z: Current position
        px, py, pz: Current momentum
        charge: Particle charge
        step_length: Step size in meters (5cm = 0.05m)
    """
    state = np.array([x, y, z, px, py, pz])
    p_mag = np.sqrt(px*px + py*py + pz*pz)
    energy = np.sqrt(p_mag*p_mag + MUON_MASS*MUON_MASS)
    
    # Calculate time step based on velocity (assuming c=1)
    # dt = ds/v where v ≈ c for ultra-relativistic particles
    #dt = step_length / (p_mag/np.sqrt(p_mag*p_mag))
    v_mag = (p_mag/energy)*c # beta = v/c = p/E
    dt = step_length / v_mag  # ds/v
    
    # RK4 integration
    k1 = derivative(state, charge, mag_field)
    k2 = derivative(state + 0.5*dt*k1, charge, mag_field)
    k3 = derivative(state + 0.5*dt*k2, charge, mag_field)
    k4 = derivative(state + dt*k3, charge, mag_field)
    
    new_state = state + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
    return tuple(new_state)

def derivative(state: np.ndarray, charge: float, mag_field):
    """
    Calculates the derivatives for position and momentum.
    """
    x, y, z, px, py, pz = state
    p_mag = np.sqrt(px*px + py*py + pz*pz)
    energy = np.sqrt(p_mag*p_mag + MUON_MASS*MUON_MASS)
    
    # Velocity components with mass (v = p/E)
    vx = (px/energy)*c
    vy = (py/energy)*c
    vz = (pz/energy)*c

    B = mag_field.get_magnetic_field(x, y, z)
    
    # Calculate force components (q[v × B])
    charge *= e_charge
    dpx = charge * (vy*B[2] - vz*B[1]) / GeV_over_c_to_SI
    dpy = charge * (vz*B[0] - vx*B[2]) / GeV_over_c_to_SI
    dpz = charge * (vx*B[1] - vy*B[0]) / GeV_over_c_to_SI
    
    return np.array([vx, vy, vz, dpx, dpy, dpz])

def track_particle(x: float, y: float, z: float,
                  px: float, py: float, pz: float,
                  charge: float, 
                  num_steps: int,
                  step_length: float = 0.05,
                  mag_field='toy'):
    """
    Tracks a particle through multiple steps with fixed spatial step size.
    
    Args:
        x, y, z: Initial position (meters)
        px, py, pz: Initial momentum (GeV/c)
        charge: Particle charge (elementary charge units)
        num_steps: Number of steps to track
        step_length: Step size in meters (default 5cm = 0.05m)
    
    Returns:
        List of (x, y, z, px, py, pz) tuples for each step
    """
    trajectory = [(x, y, z, px, py, pz)]
    
    current_x, current_y, current_z = x, y, z
    current_px, current_py, current_pz = px, py, pz
    if mag_field == 'toy':
        mag_field_generator = ToyMagneticField()
    elif mag_field == 'uniform':
        mag_field_generator = UniformMagneticField()
    else:
        mag_field_generator = CustomMagneticField(mag_field)
    for _ in range(num_steps):
        # Perform one RK4 step
        (current_x, current_y, current_z,
         current_px, current_py, current_pz) = rk4_step(
            current_x, current_y, current_z,
            current_px, current_py, current_pz,
            charge, step_length,mag_field_generator
        )
        if current_z >100:
            continue
        
        trajectory.append((
            current_x, current_y, current_z,
            current_px, current_py, current_pz
        ))
    
    return trajectory

