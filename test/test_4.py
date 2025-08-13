import numpy as np
from test_globals import *
import matplotlib.pyplot as plt
import time

# from UM_MRF import *
# import UM_MRF
# print(UM_MRF.__file__)

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from src.UM_MRF import *


if __name__ == "__main__":
    dt_dead = 1
    dt_ro = 0.1

    num_reps = 100
    RO_string = "GRE"

    PW = 2.5
    ETL = 20
    ESP = 40
    delay = 8

    sample_arr = np.array([delay + ESP + 1 + PW])
    sample_arr = np.arange(1, ETL - 10) * ESP + 2+ delay + PW
    sample_arr = np.array([delay + 10 * ESP + 1 + PW])
    crush_arr = np.arange(1, ETL) * ESP + 3.5 + delay + PW
    crush_arr = np.array([])
    print(sample_arr)
    print()
    print(crush_arr)

    #BAT = 500.0

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.000, F, lam, zvel, zpos_init, CBV, BAT, 1, 1, 10)

    # -- NON DYNAMIC TIME -- #
    ps_no_dyn = MRFSim(p)
    avg_samps_no_dyn = True

    # --- GRADIENT ECHO --- #
    # for _ in range(5):
    #     #ps_no_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=False))
    #     ps_no_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))
    #     #ps_no_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=False))
    #     ps_no_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))

    # for rep in range(num_reps):
    #     #ps_no_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=False))
    #     ps_no_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(0), dynamic_time=False))
    #     ps_no_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))
    #     #ps_no_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=False))
    #     ps_no_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(1), dynamic_time=False))
    #     ps_no_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr))


    # --- SPIN ECHO --- #
    for _ in range(5):
        #ps_no_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_no_dyn))
        #ps_no_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_no_dyn))

    for rep in range(num_reps):
        #ps_no_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(0), dynamic_time=False, avg_samples=avg_samps_no_dyn))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_no_dyn))
        #ps_no_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=False))
        ps_no_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(1), dynamic_time=False, avg_samples=avg_samps_no_dyn))
        ps_no_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=False, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_no_dyn))

    ps_no_dyn.setup()
    # ps_no_dyn.plot_B()
    # ps_no_dyn.plot_s()

    # Time it
    start_t = time.time()
    ps_no_dyn.run_all_np()
    end_t = time.time()
    print("Non-Dynamic Simulation took:  ", end_t - start_t)

    ps_no_dyn.plot_M()

    samples_no_dyn = ps_no_dyn.samples
    samples_no_dyn = samples_no_dyn[0::]
    subs_no_dyn = samples_no_dyn[0::2] - samples_no_dyn[1::2]

    plt.plot(np.abs(subs_no_dyn))
    plt.title(f"Control-Label Subtractions using the Non-Dynamic-Time simulation ({RO_string})")
    plt.xlabel("Frame Number")
    plt.ylabel("Sample Intensity [Magnetization]")
    plt.show()



    # -- DYNAMIC TIME -- #
    ps_dyn = MRFSim(p)
    avg_samps_dyn = True

    # --- GRADIENT ECHO --- #
    # for _ in range(5):
    #     #ps_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=True))
    #     ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))
    #     #ps_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=True))
    #     ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))

    # for rep in range(num_reps):
    #     #ps_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
    #     ps_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(0), dynamic_time=True))
    #     ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))
    #     #ps_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
    #     ps_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(1), dynamic_time=True))
    #     ps_dyn.add_sim(GRE(PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr))

    # --- SPIN ECHO --- #
    for _ in range(5):
        #ps_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_dyn))
        #ps_dyn.add_sim(DeadAir((num_reps - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_dyn))

    for rep in range(num_reps):
        #ps_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(0), dynamic_time=True, avg_samples=avg_samps_dyn))
        ps_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_dyn))
        #ps_dyn.add_sim(DeadAir((num_reps -rep - 1) * 100, dt_dead, dynamic_time=True))
        ps_dyn.add_sim(pCASL(rep * 100, dt_dead, control=(1), dynamic_time=True, avg_samples=avg_samps_dyn))
        ps_dyn.add_sim(FSE(500, PW, ETL, delay, ESP, dt_ro, dynamic_time=True, sample_times=sample_arr, crusher_times=crush_arr, avg_samples=avg_samps_dyn))

    ps_dyn.setup()
    # ps_dyn.plot_B()
    # ps_dyn.plot_s()

    # Time it
    start_t = time.time()
    ps_dyn.run_all_np()
    end_t = time.time()
    print("Dynamic Simulation took:      ", end_t - start_t)

    ps_dyn.plot_M()

    samples_dyn = ps_dyn.samples
    samples_dyn = samples_dyn[0::]
    subs_dyn = samples_dyn[0::2] - samples_dyn[1::2]

    plt.plot(np.abs(subs_dyn))
    plt.title(f"Control-Label Subtractions using the Dynamic-Time simulation ({RO_string})")
    plt.xlabel("Frame Number")
    plt.ylabel("Sample Intensity [Magnetization]")
    plt.show()
    

    # # Difference in the samples
    diff = subs_dyn - subs_no_dyn
    plt.plot(np.abs(diff))
    plt.title(f"Difference between the Dynamic and Non-Dynamic Simulations ({RO_string})")
    plt.xlabel("Frame Number")
    plt.ylabel("Difference in Sample Intensity [Magnetization]")
    plt.show()


    # # Error Percentage
    # perc_err = 100 * diff / subs_no_dyn
    # plt.plot(np.abs(perc_err))
    # plt.title(f"Percent error between the Dynamic and Non-Dynamic Simulations ({RO_string})")
    # plt.xlabel("Frame Number")
    # plt.ylabel("Percent Error in Sample Intensity [%]")
    # plt.show()