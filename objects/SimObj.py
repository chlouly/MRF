##########################################################################
#   This file contains the class definition for the SimObj object        #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

import numpy as np
from .simulators.np_blochsim_ljn import np_blochsim_ljn
from UM_Blochsim.blochsim import *


class SimObj:
    """
    This class represents one Block of the simulation. This codebase is structured 
    such that different different parts of a pulse sequence, like pCASL pulse trains,
    Dear Air (nothing playing), Different types of read out pulsetrains, etc. are
    represented by instances children of this class. Multiple instances of children
    of this class will be stored in an instance of MRFSim (found in './MRFSim') to be
    simulated sequentially, forming a pulse sequence.

    This is the parent class for each one of those blocks. The children classes (pCASL,
    FSE, DeadAir, etc.) can be found in dedicated files in the './sim_blocks' directory.

    Class Variables:
        absorption:         The coefficient that governs the effect of the pulse
                            on the semisolid pool. The default value is 1, however 
                            this may vary for different child classes.
        saturation:         The arbitrary constant that only serves to saturate
                            the semisolid pool. The default value is 0, however in 
                            the pCASL(SimObj) child class specifically, we use a large
                            positive value instead.
        M_start_default:    Default Starting Magnetization Vector [0 0 1 1]

    Crucial Methods:
        __init__():         Creates a new instance of the SimObj class.
        set_*():            Set's either the RF, Gradient, or arterial activation effects.
        run_np_ljn():       Runst the simulation for this block of this pulse sequence.
    """
    absorption = 1.0
    saturation = 0
    M_start_default = np.array([0.0, 0.0, 1.0, 1.0])


    def __init__(self, params, T, TR, PW, dt):
        """
        Method that instantiates the SimObj class.

        Parameters:
            params:     An instance of the Params class.
            T:          The durration of this block [ms]
            TR:         Time between pulses in the RF pulsetrains for this block [ms]
                        (Not always used)
            PW:         Pulse width [ms]
            dt:         Simulation timestep [ms]
        """
        self.T = T                                  # Simulation durration [ms]
        self.dt = dt                                # Simulation timestep [ms]
        self.ntime = int(np.ceil(T / self.dt))      # Number of time samples
        self.time_inds = np.arange(self.ntime)      # Vector of time indices
        self.time = self.time_inds * self.dt        # Vector of Timepoints [ms]

        self.B = np.zeros((self.ntime, 3))          # (ntime, 3) array of B vectors [T] (initally set to 0s)
        self.s = np.zeros((self.ntime, ))          # (ntime, ) array of arterial magnetization values (initially set to 0s)

        self.TR = TR                                
        self.PW = PW

        self.params = params


    def set_gradients(self, grads):
        """
        This incorporates the effects of the gradient fields into the total B field
        stored in an instance of the SimObj class or one of its children. It does this
        by taking the dot product of the gradient vector and the observation location 
        vector (described in this comment as r(t)) at each point in time and adding it
        to the z component of B(t):
            B_z(t) += G(t) dot r(t)
        In this method, the information about r(t) is stored in the Params object
        associated with this block.

        Parameters:
            grads:      The (ntime, 3) array of gradient vectors for each time point. 
                        [T / cm]

        Note:   This method does not get used as we are not currently considering the
                effects of gradients. This function will be left here (as well as 
                overwritten versions in child classes) to avoid future rewrites in case
                it ends up being necessary (most likely due to unforseen circumstances).
        """
        # Getting an array of z positions over time (assuming that the spin is only moving
        # through space in the z direction)
        zpos = (self.time * self.params.zvel) + self.params.zpos_init

        # Elementwise Multipication and 
        self.B[:, 2] = self.B[:, 2] + (self.params.xpos * grads[:, 0]) + (self.params.ypos * grads[:, 1]) + (zpos * grads[:, 2])


    def set_rf(self, rf):
        """
        This incorporates the effects of the RF pulses into the total B field stored 
        in an instance of the SimObj class or one of its children.

        Parameters:
            rf:     The (ntime, 3) array of B1 field vectors for each time point. 
                    [T]

        Note:   For the sake of homogenaity, B1 fields will always be described as
                3-vectors, even though the field is always in the transverse plane.
        """
        self.B = self.B + rf


    def set_s_sig(self, start_t, end_t):
        # If both queues are empty, return, nothing to do.
        if (len(start_t) == 0) or (len(end_t) == 0):
            return ([], [])
        elif start_t[0] <= self.T:
            # Pulse plays during this block
            self.s =  - (self.params.F * 2 * self.params.alpha * self.params.M0_f / self.params.lam) * np.exp(-self.params.BAT / self.params.T1_b) * ((self.time >= start_t[0]) & (self.time < end_t[0]))

        # Update start and end times
        start_t[0] = np.max((0.0, start_t[0] - self.T))
        end_t[0] = np.max((0.0, end_t[0] - self.T))

        if end_t[0] == 0:
            # If what pulse is done, we delete it
            start_t.pop(0)
            end_t.pop(0)
        
        return (start_t , end_t)
   

    def run_np_ljn(self, M_start=M_start_default):
        """
        This method runs a python implementation of the LJN Bloch simulation
        algorithm.
        """
        self.M = np_blochsim_ljn(self.B, self.s, self.params, self.dt, self.ntime, M_start, self.absorption, self.saturation, timer=True)

    

    """
    The commented sections below are methods that run simulations using an various
    Bloch Simulator implementations written in C. Currently not in use.
    """
    # def run_eul(self, M_start=M_start_default):
    #     self.M = blochsim_eul(self.B, self.params.T1_b, self.params.T2_b, self.dt)
    

    # def run_rk4(self, M_start=M_start_default):
    #     self.M = blochsim_rk4(self.B, self.params.T1_b, self.params.T2_b, self.dt,)
    

    # def run_ljn(self, M_start=M_start_default):
    #     self.M = blochsim_ljn(self.B, self.s, M_start, self.params.T1_s, self.params.T1_b, self.params.T2_b, self.dt, self.params.ks, self.params.kf, 1, self.params.F, self.params.lam)

