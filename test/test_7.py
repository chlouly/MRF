from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from test_globals import *
from objects import *

# Testing the new gradient echo block.
if __name__ == "__main__":
    dt_dead = 35
    dt_live = 0.1

    # Number of times we repeat the 4 blocks
    num_rep = 5

    # Here I am creating a Params object to give the simulation it's operating parameters
    params = Params(T1_f, T2_f, T1_s, 0.0001, 0.0001, F, lam, zvel, zpos_init, CBV, BAT, 1, 1, 90)
    sim = MRFSim(params)

    sim.add_sim(DeadAir(1000, dt_dead))
    sim.add_sim(GRE(20, 20, 100, 50, dt_live))
    #sim.add_sim(DeadAir(1000, dt_dead))

    sim.setup()

    sim.plot_B()
    sim.plot_s()

    sim.run_all_np()

    sim.plot_M()




    
