import numpy as np
from test_globals import *

from UM_MRF import MRFSim, Params, GRE
import UM_MRF
print(UM_MRF.__file__)

SCHED_PATH = "/Users/chrislouly/Documents/MICHIGAN/MRI_RESEARCH/sched_bldr"
DICT_PATH = "/Users/chrislouly/Documents/MICHIGAN/MRI_RESEARCH/dicts"


if __name__ == "__main__":
    # T1 and T2 fitting (generating dictionary)
    T1_f = np.linspace(10, 700, 10)
    T2_f = np.linspace(10, 400, 10)
    #T2_f = 1000
    flip = np.linspace(2, 20, 5)

    dyn_time = False

    PW = 2.5            # RO RF Pulse Width [ms]
    ETL = 20            # RO Echo Train Length
    ESP = 40            # RO Echo Spacing [ms]
    delay = 8           # RO Delay before echos play [ms]
    crush_off = 3.5     # RO Crusher offset (time between the end of an RF pulse and a crusher) [ms]
    samp_off = 2        # RO Sample offset (time between the end of an RF pulse and a sample) [ms]

    sample_times = (np.arange(ETL) * ESP) + delay + PW + 2
    crush_times = (np.arange(ETL) * ESP) + delay + PW + crush_off
    

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.0000, 0, lam, 0, 0, 0, BAT, 1, 1, flip)
    #p.resume(DICT_PATH + "/90007_GRE_BALL.h5")
    ps = MRFSim(p)

    ro_block = GRE(2.5, 20, 8, 40, 0.1, dynamic_time=dyn_time, crusher_times=crush_times, sample_times=sample_times, avg_samples=True)

    ps.read_sched(SCHED_PATH + "/91007", ro_block, dyn_time=dyn_time)

    ps.setup()

    ps.generate_dict(DICT_PATH + "/90007_GRE_NEW.h5")

    