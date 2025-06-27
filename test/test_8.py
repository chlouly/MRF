import numpy as np
from test_globals import *
import matplotlib.pyplot as plt

from UM_MRF import *
import UM_MRF
print(UM_MRF.__file__)

if __name__ == "__main__":
    dt_dead = 1
    dt_ro = 0.1

    num_reps = 5

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.000, F, lam, zvel, zpos_init, CBV, BAT, 1, 1, 90)
    ps = MRFSim(p)

    # for i in range(10):
    #     ps.add_sim(DeadAir((num_reps - 1) * 100, dt_dead))
    #     ps.add_sim(pCASL(2000, dt_dead, control=(i % 2)))
    #     ps.add_sim(FSE(500, 2, 20, 5, 5, dt_ro))

    for rep in range(num_reps):
        ps.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead))
        ps.add_sim(pCASL(2000, dt_dead, control=(rep % 2)))
        ps.add_sim(DeadAir(rep * 100, dt_dead))
        ps.add_sim(FSE(500, 2, 20, 5, 5, dt_ro))

    p_2 = Params(T1_f, T2_f, T1_s, 0.0000, 0.000, F, lam, zvel, zpos_init, CBV, BAT, 1, 1, 90)
    ps_2 = MRFSim(p_2)

    # for i in range(10):
    #     ps.add_sim(DeadAir((num_reps - 1) * 100, dt_dead))
    #     ps.add_sim(pCASL(2000, dt_dead, control=(i % 2)))
    #     ps.add_sim(FSE(500, 2, 20, 5, 5, dt_ro))

    for rep in range(num_reps):
        ps_2.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
        ps_2.add_sim(pCASL(2000, dt_dead, control=(rep % 2), dynamic_time=True))
        ps_2.add_sim(DeadAir(rep * 100, dt_dead, dynamic_time=True))
        ps_2.add_sim(FSE(500, 2, 20, 5, 5, dt_ro, dynamic_time=True))
        

    ps.setup()
    ps.run_all_np()

    ps_2.setup()
    ps_2.run_all_np()

    print(ps.get_times().shape)
    print(ps_2.get_times().shape)

    ps_2.setup()
    ps_2.run_all_np()

    ps.plot_B()
    ps_2.plot_B()
    
    ps.plot_s()
    ps_2.plot_s()

    ps.plot_M()
    ps_2.plot_M()

    print()
    print()

    print(ps.sims[-1].M[-1])
    print(ps_2.sims[-1].M[-1])
