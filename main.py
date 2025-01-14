"""
This program was designed to simulate inversion. 
Operating at 7T
"""
import numpy as np
from pulses import rf_gen, slice_grad_gen, time, ntime, dt
from blochsim import blochsim, x, y, z
import matplotlib.pyplot as plt

if __name__ == '__main__':
    xpos = 0.0          #[cm]
    ypos = 0.0          #[cm]
    zvel = 10**-3       #[cm/ms]
    zpos_init = -5.0    #[cm]

    zpos = (time * zvel) + zpos_init

    TR = 500 #ms
    PW = 200 #ms

    Bt = np.zeros((ntime, 3))

    Bt[:, z] = np.zeros(ntime)          # B1 field in rotating frame (z dir)
    Bt[:, x] = rf_gen(0.00005, PW, TR)     # Rf Pulse (x dir)

    Gx = np.zeros(ntime)                #
    Gy = np.zeros(ntime)                #
    Gz = slice_grad_gen(0.0003, PW, 0.1, TR)

    Bt[:, x] = Bt[:, x] + (Gx * xpos)
    Bt[:, y] = Bt[:, y] + (Gy * ypos)
    Bt[:, z] = Bt[:, z] + (Gz * zpos)

    M = blochsim(Bt)

    #plt.plot(time, rf)
    plt.plot(time * dt * 10**-3, M[:, x], label = 'x')
    plt.plot(time * dt * 10**-3, M[:, y], label = 'y')
    plt.plot(time * dt * 10**-3, M[:, z], label = 'z')
    plt.legend()
    plt.show()

