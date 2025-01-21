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


def slice_grad_gen(Gmax, pw, Gave, TR):
    # This will be one Gradient pulse
    block_len = int(np.ceil(TR / dt))
    block = np.zeros(block_len)
    p_len = int(np.ceil(pw / dt))

    # Construct the ramp up to the max amplitude
    main_ramp = np.arange(0, Gmax, max_slew * dt)
    main_ramp_len = np.size(main_ramp)

    # Construct the main gradient pulse for one repetition
    block[0:main_ramp_len] = main_ramp
    block[main_ramp_len:main_ramp_len + p_len] = Gmax
    block[main_ramp_len + p_len:2 * main_ramp_len + p_len] = main_ramp[::-1]

    # Construct the rephasing pulse
    cur_g_sum = np.sum(block) * dt
    reph_g_sum = Gave - cur_g_sum
    reph_amp = np.sqrt(np.abs(reph_g_sum * max_slew * dt))
    reph_ramp = np.arange(0, -reph_amp, -max_slew * dt)
    reph_ramp_len = np.size(reph_ramp)

    block[2 * main_ramp_len + p_len:2 * main_ramp_len + p_len + reph_ramp_len] = reph_ramp
    block[2 * main_ramp_len + p_len + reph_ramp_len:2 * main_ramp_len + p_len + 2 * reph_ramp_len] = reph_ramp[::-1]

    print(np.sum(block))

    # Create the final pulse sequence and return
    nreps = int(np.ceil(ntime / block_len))
    return np.tile(block, (nreps))[0 : ntime]
