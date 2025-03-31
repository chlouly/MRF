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

    #s = gen_s_pool_sig(plot = 0)
    s = np.zeros(ntime)

    Bt = np.zeros((ntime, n_dim))
    #Bt[:, x:z] = pcasl_rf_gen(flip, PW, TR, d_psi, psi_0, control=1)
    Bt[:, x:z] = fse_pulsetrain(PW, TR, fill=1)
    Gx = np.zeros(ntime)
    Gy = np.zeros(ntime)
    Gz = slice_grad_gen(g_max, PW, g_ave, TR)
    #Bt[0:int(np.ceil(ntime / 2)), x] = pi / (gam * 0.5)
    #Bt[:, x] = pi / (gam * 500 * 2)

    Bt[:, z] = Bt[:, z] + (Gx * xpos)
    Bt[:, z] = Bt[:, z] + (Gy * ypos)
    Bt[:, z] = Bt[:, z] + (Gz * zpos)

    print(np.shape(Bt))

    plt.plot(time[::50], Bt[::50, x], label="x")
    plt.plot(time[::50], Bt[::50, y], label="y")
    plt.xlabel("time [ms]")
    plt.ylabel("RF amplitude [T]")
    plt.title("FSE RF Pulse sequence")
    plt.legend(loc="upper right")
    plt.show()

    M = blochsim_ljn(Bt, s, np.array([0.0, 0.0, 1.0, 1.0]), T1_s, T1_b, T2_b, dt, dt, 0.0, 0.0, 1, F, lam, plot=True, dsample=1, timer=True)
    #M = blochsim_eul(Bt, T1_b, T2_b, dt, plot=True)




