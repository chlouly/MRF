"""
This program was designed to simulate inversion. 
Operating at 7T
"""
import numpy as np
from globals import *
from pulses import rf_gen, slice_grad_gen
from blochsim import blochsim
import matplotlib.pyplot as plt

if __name__ == '__main__':

    zpos = (time * zvel) + zpos_init

    Bt = np.zeros((ntime, n_dim))
    Bt[:, x:z] = rf_gen(rf_amp, PW, TR, d_psi, psi_0)

    Gx = np.zeros(ntime)
    Gy = np.zeros(ntime)
    Gz = slice_grad_gen(g_max, PW, g_ave, TR)

    Bt[:, z] = Bt[:, z] + (Gx * xpos)
    Bt[:, z] = Bt[:, z] + (Gy * ypos)
    Bt[:, z] = Bt[:, z] + (Gz * zpos)

    M = blochsim(Bt)

    #plt.plot(time_s, Gz)
    plt.plot(time_s, M[:, x], label = 'x')
    plt.plot(time_s, M[:, y], label = 'y')
    plt.plot(time_s, M[:, z], label = 'z')
    plt.legend()
    plt.show()

