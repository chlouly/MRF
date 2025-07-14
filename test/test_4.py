import numpy as np
from test_globals import *
import matplotlib.pyplot as plt
import time

from UM_MRF import *
import UM_MRF
print(UM_MRF.__file__)


if __name__ == "__main__":
    dt_dead = 1
    dt_ro = 0.1

    num_reps = 100

    PW = 2.5
    ETL = 20
    ESP = 40
    delay = 8

    sample_arr = np.array([delay + ESP + 1 + PW])
    #crush_arr = np.arange(1, ETL + 1) * ESP + 3.5 + delay + PW
    crush_arr = np.array([])

    #BAT = 500.0

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.000, F, lam, zvel, zpos_init, CBV, BAT, 1, 1, 50)

    # -- NON DYNAMIC TIME -- #
    ps_no_dyn = MRFSim(p)

    for _ in range(5):
        #ps_no_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))
        #ps_no_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))

    for rep in range(num_reps):
        #ps_no_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(0), dynamic_time=False))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))
        #ps_no_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(1), dynamic_time=False))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))

    ps_no_dyn.setup()
    #ps_no_dyn.plot_B()
    #ps_no_dyn.plot_s()
    ps_no_dyn.run_all_np()
    #ps_no_dyn.plot_M()

    samples_no_dyn = ps_no_dyn.samples
    samples_no_dyn = samples_no_dyn[10::]
    subs_no_dyn = samples_no_dyn[0::2] - samples_no_dyn[1::2]

    plt.plot(subs_no_dyn)
    plt.title("Control-Label Subtractions using the Non-Dynamic-Time simulation")
    plt.xlabel("Frame Number")
    plt.ylabel("Sample Intensity [Magnetization]")
    plt.show()



    # -- DYNAMIC TIME -- #
    ps_dyn = MRFSim(p)

    for _ in range(5):
        ps_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))
        ps_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))

    for rep in range(num_reps):
        ps_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(0), dynamic_time=True))
        ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))
        ps_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(1), dynamic_time=True))
        ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))

    ps_dyn.setup()
    #ps_dyn.plot_B()
    ps_dyn.plot_s()
    ps_dyn.run_all_np()
    #ps_dyn.plot_M()

    samples_dyn = ps_dyn.samples
    samples_dyn = samples_dyn[10::]
    subs_dyn = samples_dyn[0::2] - samples_dyn[1::2]

    plt.plot(subs_dyn)
    plt.title("Control-Label Subtractions using the Dynamic-Time simulation")
    plt.xlabel("Frame Number")
    plt.ylabel("Sample Intensity [Magnetization]")
    plt.show()
    

    # Difference in the samples
    diff = subs_dyn - subs_no_dyn
    plt.plot(diff)
    plt.title("Difference between the Dynamic and Non-Dynamic Simulations")
    plt.xlabel("Frame Number")
    plt.ylabel("Difference in Sample Intensity [Magnetization]")
    plt.show()