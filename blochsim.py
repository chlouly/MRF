import numpy as np
from pulses import dt, ntime, gam

### COORDINATE INDICES ###
x = 0
y = 1
z = 2
##########################

Minit = np.array([0, 0, 1])
Minit = np.reshape(Minit, (1, 3))

def blochsim(B, Mi=Minit):
    M = np.zeros((ntime, 3))
    M[0, :] = Mi
    print(np.cross(M[0, :], B[0, :]))

    for i in range(1, ntime):
        M[i, :] = M[i - 1, :] + gam * dt * np.cross(M[i - 1, :], B[i - 1, :])

    return M
