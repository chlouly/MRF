##########################################################################
#   This file contains the class definition for the MRFSim object        #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

import numpy as np
import matplotlib.pyplot as plt
from .dict_manip import *
from .sim_blocks import *
from .pb import create_pb, refresh_pb

M_init = np.array([0.0, 0.0, 1.0, 1.0])

class MRFSim:
    """
    This class represents a full pulse sequence. It does so by holding an array of 
    SimObj instantiations that make up the pulse sequence that we would like to simulate.

    Essentially, we break up the pulse sequence into a series of 'blocks' that are simulated
    sequentially. We use this class to wrap those blocks (SimObj instantiations) to make
    using the simulator cleaner and easier (hopefully...).

    Class Variables:
        sims:       A simple python list of SimObj instantiations. These instantiations are
                    stored in the order in which they will be simulated.
        cur_sim:    The index of the next SimObj to be simulated in the list of SimObjs.
        num_sim:    The number of SimObjs that the instantiation of this class holds. Effectively
                    the same as length(self.sims)
        M_cur:      The current magnetization vector. This represents the magnetization of the
                    SimObj that was last run, and is used as the starting magnetization for the
                    next SimObj to be run. When this class is first instantiated, M_cur is set to
                    M_cur = [0 0 1 1].

    Crucial Methods:
        __init__():     Creates a new instance of the MRFSim class.
        add_sim():      Adds a SimObj to the list of blocks to be simulated.
        setup():        Sets up each SimObj so that they are ready to be simulated.
        run_all_np():   Runs each simulation starting from the current SimObj until the end of the
                        list is reached.
    """


    def __init__(self, params):
        """
        This method creates a new and empty instantiation of MRFSim.
        """
        self.params = params
        self.cur_sim = 0
        self.num_sim = 0
        self.cur_time = 0

        self.M_cur = M_init

        self.samples = np.array([])
        self.sample_times = np.array([])

        self.sims = []


    def add_sim(self, SimObj):
        """
        This method adds a SimObj to the END of the list of blocks to be simulated.
        """
        if SimObj is None:
            # SimObj did not exist, most likely because T <= 0
            return

        self.sims.append(SimObj)
        self.num_sim += 1


    def clear(self, M_start=M_init):
        """
        Helper method that resets the simulator back to its original state.
        (Erases all SimObjs and resets all counters)
        """
        self.sims = []
        self.cur_sim = 0
        self.num_sim = 0
        self.M_cur = M_start


    def set_control(self):
        """
        Helper method that sets all pCASL blocks to be control pulsetrains
        """
        for sim in self.sims:
            sim.control = 1


    def set_label(self):
        """
        Helper method that sets all pCASL blocks to be label pulsetrains
        """
        for sim in self.sims:
            sim.control = 0


    def setup(self):
        """
        Helper method that does the following for each SimObj in the list:
            runs set_gradients() to incorporate the gradients into the sim
            runs set_rf() to incorporate the RF pulses into the sim
            runs set_s_sig() (see SimObj.set_s_sig() for details)
        """
        # Initialize an empty queue of bolus times
        time_queue = []

        T = 0.0

        for sim in self.sims:
            sim.set_gradients()
            sim.set_rf(self.params)
            if np.size(sim.sample_times) != 0:
                self.sample_times = np.append(self.sample_times, sim.sample_times + T)
            time_queue = sim.set_s_shape(time_queue, self.params.BAT)
            sim.scale_s(self.params.F, self.params.lam, self.params.alpha, self.params.M0_f, self.params.BAT, self.params.T1_b)
            T += sim.T


    def compute_s(self):
        self.params.recompute_s = False
        self.params.rescale_s = False

        time_queue = []

        for sim in self.sims:
            time_queue = sim.set_s_shape(time_queue, self.params)
            sim.scale_s(self.params.F, self.params.lam, self.params.alpha, self.params.M0_f, self.params.BAT, self.params.T1_b)

    
    def scale_s(self):
        self.params.rescale_s = False

        for sim in self.sims:
            sim.scale_s(self.params.F, self.params.lam, self.params.alpha, self.params.M0_f, self.params.BAT, self.params.T1_b)


    def modify_flips(self):
        """
        This Method runs through all simulator blocks and calls the set_flip method that changes the flip angle of the block,
        and recalculates the Effective B field. This is only the case for the gradient echo pulse, all others currently
        don't depend on a flip angle.
        
        Input Agruments:
            - flip:     Flip angle in degrees
            - phase:    RF Pulse Phase, in degrees.
        """
        for sim in self.sims:
            sim.set_flip(self.params)

    
    def read_sched(self, sched_dir: str):
        # All schedules are found in a mrf_schedule.txt file
        sched_dir = sched_dir + "/mrf_schedule.txt"

        with open(sched_dir, "r") as sched:
            for line in sched:
                # Get Values for one line
                vals = list(map(float, line.split()))

                # Add Blocks to the sim (the scheduler has times in seconds,
                # We use miliseconds).

                # First Delay
                self.add_sim(DeadAir(vals[0] * 1000, 100))

                # pCASL Section
                if int(vals[1] == -1):
                    # pCASL code of -1 means do nothing
                    self.add_sim(DeadAir(vals[2] * 1000, 300))
                else:
                    # Otherwise we do something
                    # Label == 1, Control == 0
                    self.add_sim(pCASL(vals[2] * 1000, 300, control= not int(vals[1])))
                
                # pCASL PLD
                self.add_sim(DeadAir(vals[3] * 1000, 300))

                # Prep Pulse 1
                self.add_sim(self.add_presat_sim(vals[4]))

                # Prep 1 PLD
                self.add_sim(DeadAir(vals[5] * 1000, 300))

                # Prep Pulse 2
                self.add_sim(self.add_presat_sim(vals[6]))

                # Prep 1 PLD
                self.add_sim(DeadAir(vals[7] * 1000, 300))

                # Finally We add a readout
                self.add_sim(GRE(2.5, 20, 8, 40, 2.5))


    def add_presat_sim(self, prep_pulse_code):
        code = np.int32(prep_pulse_code)
        if code == 0:     # No pulse
            return None
        elif code == 250:
            raise ValueError("Error: The prep pulse with ID 00250 has not yet been added to the MRF Simulator")
        elif code == 3200:
            raise ValueError("Error: The prep pulse with ID 03200 has not yet been added to the MRF Simulator")
        elif code == 6800:
            raise ValueError("Error: The prep pulse with ID 06800 has not yet been added to the MRF Simulator")
        elif code == 6850:  # BIR8
            return BIR8(6850 * 0.004, 0.004)
        elif code == 17268:
            raise ValueError("Error: The prep pulse with ID 17268 has not yet been added to the MRF Simulator")
        elif code == 17536:
            raise ValueError("Error: The prep pulse with ID 17536 has not yet been added to the MRF Simulator")
        elif code == 17846:
            raise ValueError("Error: The prep pulse with ID 17846 has not yet been added to the MRF Simulator")
        else:
            raise ValueError(f"Error: No prep pulse is known with ID {code:05d}")
                

    

    def run_one_np(self):
        """
        Runs the SimObj at index self.cur_sim on the list of SimObjs
        using the LJN simulator written in python.
        """
        # Run the current SimObj
        #self.sims[self.cur_sim].run_np_ljn(self.params, self.M_cur)
        self.sims[self.cur_sim].run_ljn(self.params, self.M_cur)

        # Take the LAST vector from the simulated magnetization
        # and use it as the starting magnetization for the next sim (M_cur)
        self.M_cur = self.sims[self.cur_sim].M[-1, :]

        # Get the samples from the current simulation and ass them to the array
        # of all samples, as well as the time in which they were taken w.r. to
        # the pulse sequence
        if np.size(self.sims[self.cur_sim].sample_times) != 0:
            self.samples = np.append(self.samples, self.sims[self.cur_sim].sample(self.params.CBV))
            #self.sample_times = np.append(self.sample_times, self.sims[self.cur_sim].sample_times + self.cur_time)

        # Increment the current time
        self.cur_time += self.sims[self.cur_sim].T

        # Increment the index for the current simulation
        self.cur_sim += 1


    def run_all_np(self):
        """
        Runs all simulations in succession starting from the current
        SimObj.
        """
        while self.cur_sim < self.num_sim:
            self.run_one_np()


    # FOR NEXT COMMIT
    def generate_dict(self, dict_filename):
        # Initialize params so that we can iterate over it
        iter(self.params)

        # Initialize the dictionary file
        init_dict(dict_filename, self.params, np.size(self.sample_times))

        # Create Progress Bar
        create_pb()

        # Do the actual looping now
        try:
            while True:
                # Modify s(t) if needed
                if self.params.recompute_s:
                    self.compute_s()
                elif self.params.rescale_s:
                    self.scale_s()

                if self.params.recompute_B:
                    self.modify_flips()

                # Run simulations for the entire pulse sequence
                self.run_all_np()

                #Store samples
                store_entry(dict_filename, self.params.get_cur_idx(), self.samples)

                # Soft reset to prepare for the next run
                self.soft_reset()

                # Refresh the progress bar
                refresh_pb(self.params.get_comp_perc())

                # Move on to the next set of parameters
                next(self.params)

        except StopIteration:
            print("Dictionary Generation Complete!!")


    def soft_reset(self):
        self.cur_sim = 0
        self.cur_time = 0
        self.M_cur = M_init
        self.samples = np.array([])
        #self.sample_times = np.array([])


    def get_times(self):
        """
        Helper function that returns an array of time points accross all
        SimObj blocks. The time array from each SimObj are combined in order
        to create an array of all sample times throughout the entire simulation.

        This is mainly used for plotting purposes. Once this array is calculated, 
        it is stored in case it is needed later.
        """

        if hasattr(self, "time_vec"):
            # We have already calculated this, nothing to do.
            return self.time_vec
        
        cur_time = 0.0
        self.time_vec = np.empty(0)

        for sim in self.sims:
            self.time_vec = np.append(self.time_vec, sim.time + cur_time)
            cur_time += sim.T
    

    def get_M(self):
        """
        Helper method that concatenates the magnetization from each SimObj block
        sequentially.

        This is mainly used for plotting purposes. Once this array is calculated, 
        it is stored in case it is needed later.
        """
        M_out = np.empty((0, 4))

        for sim in self.sims:
            M_out = np.append(M_out, sim.M, axis=0)

        return M_out
    

    def get_B(self):
        """
        Helper method that concatenates the effective B field array from each SimObj
        block sequentially.

        This is mainly used for plotting purposes. Once this array is calculated, 
        it is stored in case it is needed later.
        """
        B_out = np.empty((0, 3))

        for sim in self.sims:
            B_out= np.append(B_out, sim.B, axis=0)

        return B_out


    def get_s(self):
        """
        Helper method that concatenates the arterial magnetization array from each 
        SimObj block sequentially.

        This is mainly used for plotting purposes. Once this array is calculated, 
        it is stored in case it is needed later.
        """
        s_out = np.empty((0,))

        for sim in self.sims:
            s_out= np.append(s_out, sim.s, axis=0)

        return s_out


    def plot_B(self, dsample=1, ylim=[]):
        """
        Helper function that plots the effective B field for the entire simulation,
        allowing users to display the pulse sequence that they are working with.
        """

        # Get the things to plot
        self.get_times()
        B = self.get_B()

        plt.plot(self.time_vec[::dsample ], B[::dsample, 0], label = 'x')
        plt.plot(self.time_vec[::dsample ], B[::dsample, 1], label = 'y')
        plt.plot(self.time_vec[::dsample ], B[::dsample, 2], label = 'z')
        if not ylim == []:
            plt.ylim(ylim)
        plt.xlabel("Time [ms]")
        plt.ylabel("B (T)")
        plt.title("Effective B field")
        plt.legend()
        plt.show()


    def plot_M(self, dsample=1, ylim=[]):
        """
        Helper function that plots the simulated magnetization for the entire simulation,
        allowing users to display results of the simulation.
        """

        # Get the things to plot
        self.get_times()
        M = self.get_M()

        plt.plot(self.time_vec[::dsample ], M[::dsample, 0], label = 'x tissue')
        plt.plot(self.time_vec[::dsample ], M[::dsample, 1], label = 'y tissue')
        plt.plot(self.time_vec[::dsample ], M[::dsample, 2], label = 'z tissue')
        plt.plot(self.time_vec[::dsample ], M[::dsample, 3], label = 'semisolid')
        if not ylim == []:
            plt.ylim(ylim)
        plt.xlabel("Time [ms]")
        plt.ylabel("Magnetization ")
        plt.title("Simulated Magnetization [T/volume]")
        plt.legend()
        plt.show()


    def plot_s(self, dsample=1, ylim=[]):
        """
        Helper function that plots arterial magnetization for the entire simulation.
        """

        # Get the things to plot
        self.get_times()
        s = self.get_s()

        plt.plot(self.time_vec[::dsample ], s[::dsample], label = 's(t)')
        if not ylim == []:
            plt.ylim(ylim)
        plt.xlabel("Time [ms]")
        plt.ylabel("Arterial magnetization (z component)")
        plt.title("s(t) signal")
        plt.legend()
        plt.show()

    def plot_samples(self, ylim=[]):
        """
        Helper function that plots the samples of the pulse sequence
        """

        if np.size(self.samples) == 0:
            # We did not collect any samples so there isn't anything
            # to do
            return
        
        plt.plot(self.sample_times, self.samples)
        if not ylim == []:
            plt.ylim(ylim)
        plt.xlabel("Sample Time [ms]")
        plt.ylabel("Sample Intensity")
        plt.title("Samples")
        plt.show()

