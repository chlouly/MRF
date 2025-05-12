from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from objects import *
import matplotlib.pyplot as plt


if __name__ == "__main__":
    dt_dead = 40
    dt_ro = 0.1

    num_reps = 100

    ps = MRFSim()
    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.000, F, lam, zvel, zpos_init, CBV, BAT, 1, 1)

    for _ in range(10):
        #ps.add_sim(DeadAir(p, (num_reps - 1) * 100, dt_dead))
        ps.add_sim(FSE(p, 500, 5, 2, dt_ro))
        #ps.add_sim(DeadAir(p, (num_reps - 1) * 100, dt_dead))
        ps.add_sim(FSE(p, 500, 5, 2, dt_ro))

    for rep in range(num_reps):
        #ps.add_sim(DeadAir(p, (num_reps -rep - 1) * 100, dt_dead))
        ps.add_sim(pCASL(p, rep * 100, dt_dead, control=(0)))
        ps.add_sim(FSE(p, 500, 5, 2, dt_ro))
        #ps.add_sim(DeadAir(p, (num_reps -rep - 1) * 100, dt_dead))
        ps.add_sim(pCASL(p, rep * 100, dt_dead, control=(1)))
        ps.add_sim(FSE(p, 500, 5, 2, dt_ro))
        

    ps.setup()
    ps.run_all_np()
    #ps.plot_samples()

    samples = ps.samples
    samples = samples[20::]
    subs = samples[0::2] - samples[1::2]

    # M_art = np.array([])
    # for sim in ps.sims:
    #     M_art = np.append(M_art, sim.M_art)

    plt.plot(subs)
    plt.show()
    