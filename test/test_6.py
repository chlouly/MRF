from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from objects import *


if __name__ == "__main__":
    # T1 and T2 fitting
    T1_f = np.linspace(0, 1000, 100)
    T2_f = np.linspace(0, 2000, 100)
    

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.0000, 0, lam, 0, 0, CBV, BAT, 1, 1)
    ps = MRFSim(p)

    # Read ps schedule file to initialize

    ps.generate_dict("90002.h5")

    