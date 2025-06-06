from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from objects import *

SCHED_PATH = "/home/clouly/DEV/sched_bldr"
DICT_PATH = "/home/clouly/DEV/dicts"


if __name__ == "__main__":
    # T1 and T2 fitting (generating dictionary)
    T1_f = np.linspace(10, 1000, 100)
    T2_f = np.linspace(10, 2000, 100)
    

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.0000, 0, lam, 0, 0, CBV, BAT, 1, 1)
    ps = MRFSim(p)

    ps.read_sched(SCHED_PATH + "/90002")

    ps.setup()

    ps.generate_dict(DICT_PATH + "/90002.h5")

    