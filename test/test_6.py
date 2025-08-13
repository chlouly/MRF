# import numpy as np
# from test_globals import *

# from UM_MRF import MRFSim, Params, GRE
# import UM_MRF
# print(UM_MRF.__file__)

from pathlib import Path
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from src.UM_MRF import *
import time

SCHED_PATH = "/home/clouly/DEV/sched_bldr"
DICT_PATH = "/home/clouly/DEV/dicts"


if __name__ == "__main__":
    # T1 and T2 fitting (generating dictionary)
    T1_f = np.linspace(400, 2000, 100)
    T2_f = np.linspace(1, 200, 100)
    #T2_f = 20
    flip = np.linspace(10, 30, 10)
    #flip = 50 * 0.67

    dyn_time = True

    PW = 2.5            # RO RF Pulse Width [ms]
    ETL = 20            # RO Echo Train Length
    ESP = 40            # RO Echo Spacing [ms]
    delay = 8           # RO Delay before echos play [ms]
    crush_off = 3.5     # RO Crusher offset (time between the end of an RF pulse and a crusher) [ms]
    samp_off = 2        # RO Sample offset (time between the end of an RF pulse and a sample) [ms]

    sample_times = (np.arange(ETL) * ESP) + delay + PW + 2
    #crush_times = (np.arange(ETL) * ESP) + delay - crush_off
    crush_times = np.array([])

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.0000, 0, lam, 0, 0, 0, BAT, 1, 1, flip)
    #p.resume(DICT_PATH + "/90007_GRE_BALL.h5")
    ps = MRFSim(p)

    ro_block = GRE(2.5, 20, 8, 40, 1, 15, dynamic_time=dyn_time, crusher_times=crush_times, sample_times=sample_times, avg_samples=True)
    
    # Create Prep Pulse Block
    rho = np.loadtxt("../pulses/06850/rho.txt", dtype=np.float64)[:, 0]
    #theta = np.loadtxt("../pulses/06850/theta.txt", dtype=np.float64)[:, 0] * np.pi / (2**15)
    theta = np.zeros((len(rho),))
    print(np.shape(theta))

    rho *= (234 * (10 ** -4)) / np.max(rho)

    B_prep1 = np.zeros((len(rho), 3))
    B_prep1[:, 0] += rho * np.cos(theta)
    B_prep1[:, 1] += rho * np.sin(theta)

    s_prep1 = np.zeros((len(rho), 1))

    prep1_block = Custom(B_prep1, s_prep1, 0.004, dynamic_time = dyn_time, crusher_times=np.array([6450 * 0.004]))



    ps.read_sched(SCHED_PATH + "/90008", ro_block, prep1_block, None, dyn_time=dyn_time)

    ps.concatenate_deadtimes()

    ps.setup()

    start_t = time.time()
    ps.generate_dict(DICT_PATH + "/90008_REAL_BIR8_BIGRHO.h5")
    end_t = time.time()

    print(f"Dictionary Generation Finished In: {end_t - start_t} seconds")
    






    # #1
    # p_1 = Params(T1_f[50], T2_f[10], T1_s, 0.0000, 0.0000, 0, lam, 0, 0, 0, BAT, 1, 1, flip[10])
    # #p.resume(DICT_PATH + "/90007_GRE_BALL.h5")
    # ps_1 = MRFSim(p_1)

    # ro_block = GRE(2.5, 20, 8, 40, 0.1, dynamic_time=dyn_time, crusher_times=crush_times, sample_times=sample_times, avg_samples=True)

    # ps_1.read_sched(SCHED_PATH + "/90008", ro_block, dyn_time=dyn_time)

    # ps_1.setup()

    # # ps.generate_dict(DICT_PATH + "/90008_FLOW_corr_crush.h5")

    # ps_1.run_all_np()


    # #2
    # p_2 = Params(T1_f[50], T2_f[10], T1_s, 0.0000, 0.0000, 0, lam, 0, 0, 0, BAT, 1, 1, flip[29])
    # #p.resume(DICT_PATH + "/90007_GRE_BALL.h5")
    # ps_2 = MRFSim(p_2)

    # ro_block = GRE(2.5, 20, 8, 40, 0.1, dynamic_time=dyn_time, crusher_times=crush_times, sample_times=sample_times, avg_samples=True)

    # ps_2.read_sched(SCHED_PATH + "/90008", ro_block, dyn_time=dyn_time)

    # ps_2.setup()

    # # ps.generate_dict(DICT_PATH + "/90008_FLOW_corr_crush.h5")

    # ps_2.run_all_np()

    # print(T1_f[50], T2_f[10], flip[3])



    # samps_1 = ps_1.samples
    # samps_2 = ps_2.samples

    # plt.plot(np.abs(samps_1), label=f"flip={flip[10]}")
    # plt.plot(np.abs(samps_2), label=f"flip={flip[29]}")
    # plt.legend()
    # plt.show()
    