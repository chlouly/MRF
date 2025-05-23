from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from objects import *

"""
This is the most basic test case. It is intended to test the validity of the actual simulations
not necessarilly the infrastructure around them. Run this test case and visually inspect the
resulting figures. This consists of 3 separate simulations, which are as follows:

    1)  An inversion pulse in x
    2)  An inversion pulse in y
    3)  Dead time to observe the decay of the signal.

Good luck!
"""

if __name__ == "__main__":
    # Test to visually confirm that simulation outputs are correct
    dt = 10**-1
    T = 1000
    ntime = int(np.ceil(T / dt))

    time = np.arange(ntime)
    s_empty = np.zeros(ntime)


    #  -  Test 1, 180x pulse  -  #
    p1 = Params(1000000, 1000000, 1000000, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1)
    Test_sim = MRFSim(p1)
    B1 = np.zeros((ntime, 3))
    B1[:, 0] = pi / (gam * T)
    Test_sim.add_sim(Custom(B1, s_empty, dt))

    Test_sim.plot_B()
    Test_sim.plot_s()

    Test_sim.setup()
    Test_sim.run_all_np()
    Test_sim.plot_M()
    Test_sim.clear()


    #  -  Test 2, 180y pulse  -  #
    p2 = p1
    Test_sim = MRFSim(p2)
    B2 = np.zeros((ntime, 3))
    B2[:, 1] = pi / (gam * T)
    Test_sim.add_sim(Custom(B2, s_empty, dt))

    Test_sim.plot_B()
    Test_sim.plot_s()

    Test_sim.setup()
    Test_sim.run_all_np()
    Test_sim.plot_M()
    Test_sim.clear(M_start=np.array([1.0, 0.0, -1.0, -1.0]))


    #  -  Test 3, Decay  -  #
    p3 = Params(100, 300, 1000, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1)
    Test_sim = MRFSim(p3)
    B3 = np.zeros((ntime, 3))
    Test_sim.add_sim(Custom(B3, s_empty, dt))

    Test_sim.plot_B()
    Test_sim.plot_s()

    Test_sim.setup()
    Test_sim.run_all_np()
    Test_sim.plot_M()
    Test_sim.clear()
