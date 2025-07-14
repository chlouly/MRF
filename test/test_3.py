import numpy as np
from test_globals import *
import matplotlib.pyplot as plt

from UM_MRF import *
import UM_MRF
print(UM_MRF.__file__)


if __name__ == "__main__":
    dt_dead = 40
    dt_ro = 0.1

    ETL = 20
    ESP = 5
    PW = 2
    delay = 5

    num_reps = 100

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.000, F, lam, zvel, zpos_init, CBV, 1500, 1, 1, 90)
    ps = MRFSim(p)

    for i in range(10):
        ps.add_sim(DeadAir((num_reps - 1) * 100, dt_dead))
        ps.add_sim(pCASL(2000, dt_dead, control=(i % 2)))
        ps.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, sample_times=np.array([ESP + 1])))

    for rep in range(num_reps):
        ps.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead))
        ps.add_sim(pCASL(2000, dt_dead, control=(rep % 2)))
        ps.add_sim(DeadAir(rep * 100, dt_dead))
        ps.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, sample_times=np.array([ESP + 1])))
        

    ps.setup()
    ps.run_all_np()
    #ps.plot_samples()

    samples = ps.samples
    samples = samples[::]
    subs = samples[0::2] - samples[1::2]

    plt.plot(subs)
    plt.show()
    