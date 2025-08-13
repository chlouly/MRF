import numpy as np
from test_globals import *
import matplotlib.pyplot as plt

from UM_MRF import *
import UM_MRF
print(UM_MRF.__file__)
from copy import deepcopy



if __name__ == '__main__':
    """
    This is currently just to test one set of pulses. This will all be replaced when the infrastructure to generate
    full dictionaries (basically just doing this many many times over) is complete.

    Here we are simulating the following sequence:

        | Dead Air | pCASL | Dead Air | FSE RO | 

    Once with a label, and once with a control. We then subtract the two to see if we get the desired output
    """
    dt_dead = 35
    dt_live = 0.1

    # Number of times we repeat the 4 blocks
    num_rep = 5

    # Here I am creating a Params object to give the simulation it's operating parameters
    params = Params(T1_f, T2_f, T1_s, 0.0001, 0.0001, F, lam, zvel, zpos_init, CBV, BAT, 1, 1, 90)

    # Here, I create an instance of one simulator object (this one is the labeling case)
    sim_l = MRFSim(params)

    # Now I add all 4 blocks of the sequence
    for _ in range(num_rep):
        sim_l.add_sim(DeadAir(2000, dt_dead))               # Add dead air
        sim_l.add_sim(pCASL(2000, dt_live))         # Add the first pCASL sequence
        sim_l.add_sim(DeadAir(5000, dt_dead))               # Add the second dead air
        sim_l.add_sim(FSE(500, 3, 20, 10, 40, dt_live))             # Add the readout section

    # I deepcopy to create an actual copy of the labeling simulator
    # This will be a separate object with all of the same elements as the first one
    # This simulator object will be the controll simulator
    sim_c = deepcopy(sim_l)         

    # Now I set each of the simulators to eather label or control
    sim_l.set_label()               # Set to label
    sim_c.set_control()             # Set to control

    # These method calls "setup" the two simulators by...
    #   Calculating and setting their RF pulses
    #   Calculating and setting their gradients (not currently in use...)
    #   Calculating and setting their s(t) functions
    sim_l.setup()
    sim_c.setup()

    # General plorring to validate that all of the state arrays (B(t), s(t)) make sense
    sim_l.plot_B(dsample=10)
    sim_l.plot_s(dsample=10)
    sim_c.plot_B(dsample=10)
    sim_c.plot_s(dsample=10)

    # Run all blocks (using the simulator written in python with numpy)
    sim_l.run_all_np()
    sim_c.run_all_np()

    # Plot the simulated magnetizations
    #sim_l.plot_M(dsample=10)
    #sim_c.plot_M(dsample=10)

    sim_l.plot_samples()
    sim_c.plot_samples()

    # Subtract the label and control magnetizations
    M_sub = sim_l.get_M() - sim_c.get_M()
    time_vec = sim_l.get_times()
    

    # Plot the sibtraction
    ylim = []
    dsample = 10
    
    plt.plot(time_vec[::dsample ], M_sub[::dsample, 0], label = 'x tissue')
    plt.plot(time_vec[::dsample ], M_sub[::dsample, 1], label = 'y tissue')
    plt.plot(time_vec[::dsample ], M_sub[::dsample, 2], label = 'z tissue')
    plt.plot(time_vec[::dsample ], M_sub[::dsample, 3], label = 'semisolid')
    if not ylim == []:
        plt.ylim(ylim)
    plt.xlabel("Time [ms]")
    plt.ylabel("Magnetization ")
    plt.title("Simulated Magnetization [T/volume]")
    plt.legend()
    plt.show()