import numpy as np

T = 10000                   # Simulation durration [ms]
dt = 10**-2                      # Simulation timestep [ms]
ntime = int(np.ceil(T / dt))   # Number of time samples
gambar = 42570              # Gyromagnetic coefficient [kHz/T]
gam = gambar * 2 * 3.14159       # Gamma [kRad/sT]

time = np.arange(ntime) * dt # Time vector [ms]

def rf_gen(amp, pw, TR, control=False):
    # This will be one RF pulse
    block = np.zeros(int(np.ceil(TR / dt)))
    p_len = int(np.ceil(pw / dt))

    # This is our pulse shape for one pulse. We are using hanning
    block[0:p_len] = amp * np.hanning(p_len)

    # Combine them into a full pulse sequence:
    if control:
        # Alternating + and - RF pulses for the Control
        superblock = np.concatenate([block, -1 * block])
    else:
        # Only + signals are used for the Labeling pulses
        superblock = np.concatenate([block, block])

    # Create the final pulse sequence and return
    nreps = int(np.ceil(ntime / (2 * p_len)))
    return np.tile(superblock, (nreps))[0 : ntime]


def slice_grad_gen(amp, pw, Gave, TR):
    # This will be one Gradient pulse
    block = np.zeros(int(np.ceil(TR / dt)))
    p_len = int(np.ceil(pw / dt))

    rephase_fact = 4    # How many times longer the main puse is than the rephasing pulse
    p_rephase_len = int(np.minimum(np.ceil(p_len / rephase_fact), np.floor(TR / dt)))

    block[0:p_len] = amp
    block[p_len:p_len + p_rephase_len] = Gave - amp * rephase_fact

    # Create the final pulse sequence and return
    nreps = int(np.ceil(ntime / p_len))
    return np.tile(block, (nreps))[0 : ntime]
