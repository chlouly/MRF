"""
This program was designed to simulate inversion. 
Operating at 7T
"""
import numpy as np
from globals import *
from pcasl import *
from fse import *
from UM_Blochsim.blochsim import *
import matplotlib.pyplot as plt

if __name__ == '__main__':
    zpos = (time * zvel) + zpos_init

    Bt = np.zeros((ntime, n_dim))
    # Bt[:, x:z] = rf_gen(rf_amp, PW, TR, d_psi, psi_0, control=0)
    Bt[:, x:z] = fse_pulsetrain(PW, TR, fill=True)
    #Bt[:, x] = Bt[:, x] + 0.0001

    Gx = np.zeros(ntime)
    Gy = np.zeros(ntime)
    Gz = slice_grad_gen(g_max, PW, g_ave, TR)

    Bt[:, z] = Bt[:, z] + (Gx * xpos)
    Bt[:, z] = Bt[:, z] + (Gy * ypos)
    Bt[:, z] = Bt[:, z] + (Gz * zpos)

    plt.plot(time[::50], Bt[::50, x], label="x")
    plt.plot(time[::50], Bt[::50, y], label="y")
    plt.xlabel("time [ms]")
    plt.ylabel("RF amplitude [T]")
    plt.title("FSE RF Pulse sequence")
    plt.legend(loc="upper right")
    plt.show()

    M = blochsim_rk4(Bt, T1, T2, dt, plot=True, dsample=50, timer=True)




