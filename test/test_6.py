import numpy as np
from test_globals import *

from UM_MRF import MRFSim, Params
import UM_MRF
print(UM_MRF.__file__)

SCHED_PATH = "/Users/chrislouly/Documents/MICHIGAN/MRI_RESEARCH/sched_bldr"
DICT_PATH = "/Users/chrislouly/Documents/MICHIGAN/MRI_RESEARCH/dicts"


if __name__ == "__main__":
    # T1 and T2 fitting (generating dictionary)
    T1_f = np.linspace(10, 700, 10)
    T2_f = np.linspace(10, 400, 10)
    #T2_f = 1000
    flip = np.linspace(1, 20, 1)
    

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.0000, 0, lam, 0, 0, 0, BAT, 1, 1, flip)
    #p.resume(DICT_PATH + "/90007_GRE_BALL.h5")
    ps = MRFSim(p)

    ps.read_sched(SCHED_PATH + "/91007")

    ps.setup()

    ps.generate_dict(DICT_PATH + "/90007_GRE_BALL_better_T2.h5")

    