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


    def __new__(cls, params, T, *args, **kwargs):
        # Basically just input validation, not important
        if T <= 0:
            # If we have 0 or negative time, nothing happens so we do not create an instance
            return None
        return super().__new__(cls)


    def __init__(self, params, T, TR, PW, dt, sample_times=np.array([])):
        """
        Method that instantiates the SimObj class.

        Parameters:
            params:         An instance of the Params class.
            T:              The durration of this block [ms]
            TR:             Time between pulses in the RF pulsetrains for this block [ms]
                            (Not always used)
            PW:             Pulse width [ms]
            dt:             Simulation timestep [ms]
            sample_times:   An array of times at which you would like to sample the signal.
                            These times must be within the time range of the simulation
                            block [0, T) ms. They will be rounded down to the nearest 
                            time value that is simulated.
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

        # Input validation: Make sure the sample points are within the block
        if np.any((sample_times >= self.T) | (sample_times < 0.0)):
            raise ValueError("ERROR: Sample times must be within the length of the block [0.0, ", self.T, ") ms.")

        self.sample_times = sample_times
        self.M_art = np.array([])   # Array of samples of the arterial magnetization z component M_art(t)

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


    def set_s_sig(self, time_queue):
        """
        This method sets the s(t) function for each block. The MRFSim object calls
        this method for each block. If a block is a pCASL labeling block, then a start
        and end time are added to the time_queues (whos state is held in
        the MRFSim.setup() method). These start and end times denote the beginning
        and end of each rect function in s(t).

        Input Parameters:
            time_queue:     A queue of start times and end times for incoming boluses. The times
                            in [ms] are w.r. to the starting time of the block, i.e. if a block 
                            begins at time 2000 ms, a start time of 1000 ms would mean 1000 ms AFTER the
                            beginning of the block (3000 ms). Similarly if the end time was 2000 ms,
                            then the bolus would pass at 4000 ms with respect to the pulse sequence.
                            This queue is a list of the following tuples:
                                (Bolus Arrival [ms], Bolus Departure [ms])
                            Where each entry refers to one bolus of labeled blood.
        
        Output Values:
            time_queue:     An updated version of the input queue. All times in the queue are shifted
                            to be w.r. to the NEXT block in the sequence, and if a bolus has already
                            passed by that point, it is removed from the queue. This state is passed 
                            off to the calling function (MRFSim.setup()), so that it may be used as 
                            input when this method is called on the next block in the pulse sequence 
                            (MRFSim.sims).
        
        NOTE: If set_s_sig() is called on an instantiation of the pCASL(SimObj) class, the control
              flow of this program is handed off to its definition of this method instead of directly
              to this definition. The pCASL(SimObj).set_s_sig() method essentially adds a start time
              and an end time to the end of the queue ONLY if it is a labeling block, not
              control. After it does that, it calls the SimObj.set_s_sig() and passes it's modified
              queue as input, essentially making it so that upon seeing a labeling sequence, the
              code knows that a new bolus is coming.
        """
        # If the queue is empty, return, nothing to do.
        if len(time_queue) == 0:
            return ([])
        elif time_queue[0][0] < self.T:
            # Pulse plays during this block iff the start time of the pulse is less
            # than the durration of the block
            self.s =  - (self.params.F * 2 * self.params.alpha * self.params.M0_f / self.params.lam) * \
                np.exp(-self.params.BAT / self.params.T1_b) * ((self.time >= time_queue[0][0]) & (self.time < time_queue[0][1]))
            
            # We also sample the arterial magnetization M_art
            self.M_art = 1 - (2 * self.params.alpha * np.exp(-(self.sample_times - time_queue[0][0]) / self.params.T1_b)) \
                * ((self.sample_times >= time_queue[0][0]) & (self.sample_times < time_queue[0][1]))

        # Update start and end times
        # We now want to make the start and end times w.r. to the 
        # beginning of the next block, so we subtract the length of
        # this block.
        #
        # If a start time is < 0, we know that part of
        # the bolus appears in this block.
        # If an end time is <= 0, we know that bolus has passed, so
        # we can take it off the queue
        #
        # The following line of code performs the logic described above
        # and returns the queue:
        return [(t[0] - self.T, t[1] - self.T) for t in time_queue if t[1] > self.T]
   

    def run_np_ljn(self, M_start=M_start_default):
        """
        This method runs a python implementation of the LJN Bloch simulation
        algorithm.
        """
        self.M = np_blochsim_ljn(self.B, self.s, self.params, self.dt, self.ntime, M_start, self.absorption, self.saturation, timer=False)

    
    def sample(self):
        """
        This function obtains samples from the simulated magnetization and returns them. This
        method will be called from the MRFSim object where the samples will be stored and added
        to a dictionary. We use the following to calculate samples:

            Samples = (1 - CBV) * |M_xy(t_smaple)|_l2 + CBV * M_art(t_sample) * sin(beta)

        Where beta is the flip angle from the RO sequence. Since we are simulating the effects of 
        the readout on the tissue magnetization (first 3 of the 4 components of M(t) or self.M in 
        the SimObj instance),we do not need to use sin(beta) to emulate the readout, we just take the
        l2 norm of the transverse components. However we do not simulate the readout on the arterial
        magnetization (not the same as s(t)), we use the sin(beta) to capture the RO effects on the blood.

        NOTE: We are using beta = 90  deg here for now (sin(90) = 1, so it isn't written in the code), 
        since our readouts at the moment start with a 90y pulse, however a more general approach may 
        be used in the future if needed.

        Output:
            An (n, ) numpy array of samples from this block, where n is the number of sample times.
        """

        # Get an array of time indices for each sample
        # TODO: If we have more than one sample per block, worry about the order (ascending time)
        self.sample_inds = np.int32(np.floor(self.sample_times / self.dt))

        # return np.linalg.norm(self.M[self.sample_inds, 0:2], axis=1) * (1 - self.params.CBV) \
        #         + self.M_art * self.params.CBV


        return np.linalg.norm(self.M[self.sample_inds, 0:2], axis=1) * (1 - self.params.CBV) \
                + self.s[self.sample_inds] * self.params.CBV


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

