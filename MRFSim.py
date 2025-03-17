import numpy as np
from UM_Blochsim.blochsim import *

class MRFSim:
    def __init__(self):
        self.cur_sim = 0
        self.num_sim = 0

        self.sims = []

    def add_sim(self, SimObj):
        self.sims.append(SimObj)
        self.num_sims = self.num_sims + 1

    def run_one_rk4(self):
        if (self.num_sim <= 0) | (self.cur_sim > self.num_sim):
            # Nothing to simulate
            return
        
        M = blochsim_rk4(self.B[self.cur_B], )
        
        self.cur_B = self.cur_B + 1
        self.cur_M = self.cur_M + 1
