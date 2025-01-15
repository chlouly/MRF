import numpy as np
from globals import *
import matplotlib.pyplot as plt


def rf_gen(amp, pw, TR, d_psi, psi_0, control=False):
    # This will be one RF pulse
    block = np.zeros(int(np.ceil(TR / dt)))
    p_len = int(np.ceil(pw / dt))

    # This is our pulse shape for one pulse. We are using hanning
    block[0:p_len] = amp * np.hanning(p_len)
    b_len = np.size(block)

    # For control, we add phi to the phase 
    if control:
        d_psi = d_psi + pi
    
    # Create the final pulse sequence and return
    n_reps = int(np.ceil(ntime / b_len))
    B_out = np.zeros((n_reps * b_len, 2))

    for i in range(n_reps):
        B_out[i * b_len : (i + 1) * b_len, x] = block * np.cos(psi_0 + i * d_psi)
        B_out[i * b_len : (i + 1) * b_len, y] = block * np.sin(psi_0 + i * d_psi)

    return B_out[0:ntime, :]


def slice_grad_gen(amp, pw, Gave, TR):
    # This will be one Gradient pulse
    block = np.zeros(int(np.ceil(TR / dt)))
    p_len = int(np.ceil(pw / dt))

    p_rephase_len = int(np.minimum(np.ceil(p_len / rephase_fact), np.floor(TR / dt)))

    block[0:p_len] = amp
    block[p_len:p_len + p_rephase_len] = Gave - amp * rephase_fact

    # Create the final pulse sequence and return
    nreps = int(np.ceil(ntime / p_len))
    return np.tile(block, (nreps))[0 : ntime]
