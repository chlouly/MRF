##########################################################################
#   This file contains the class definition for the Params object        #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################

import numpy as np
from numbers import Number

class Params:
    """
    This class holds all simulation parameters that are not related layout and timing
    of the pulse sequence. 
    """


    def __init__(self, T1_f, T2_f, T1_s, ks, kf, F, lam, zvel, zpos_init, CBV, BAT, M0_f, M0_s, alpha = 0.86):
        """
        This method initializes an instance of the Params object. 

        Iterable Input Parameters:
            T1_f:       Free water (tissue) T1 [ms]
            T2_f:       Free water (tissue) T2 [ms]
            T1_s:       Semisolid pool (spins bound by macro-molecules) T1 [ms]
            ks:         Magnetization transfer rate (Semisolid -> Free Water) [ms]
            kf:         Magnetization transfer rate (Free Water -> Semisolid) [ms]
            F:          Flow Rate [mL / g / ms] ~= [1 / ms]
            CBV:        Cerebral Blood Volume [unitless]
            BAT:        Bollus Arrival Time [ms]
            alpha:      Labelling efficiency [unitless]

        Non-Iterable Input Parameters:
            lam:        Blood Brain Barrier Coeff [mL / g]
            zvel:       Observation Spin z velocity [cm / ms]
            zpos_init:  Observation spin initial z position [cm]
            M0_f:       Longitudinal Equilibrium Magnetization of the free water pool
            M0_s:       Longitudinal Equilibrium Magnetization of the semisolid pool
        """
        #TODO: update the input params comment
        # ARRAYS OF VALUES
        self.T1_f_vals = arr_or_num(T1_f)
        self.T2_f_vals = arr_or_num(T2_f)
        self.T1_s_vals = arr_or_num(T1_s)
        self.ks_vals = arr_or_num(ks)
        self.kf_vals = arr_or_num(kf)
        self.F_vals = arr_or_num(F)
        self.CBV_vals = arr_or_num(CBV)
        self.BAT_vals = arr_or_num(BAT)
        self.alpha_vals = arr_or_num(alpha)

        # Defining simulation constants
        self.T1_b = 1600                        # Blood T1 (set to the typical 1600 ms)   
        self.lam = lam
        self.M0_f = M0_f
        self.M0_s = M0_s
        self.f = M0_s / M0_f

        self.xpos = 0                           #[cm]
        self.ypos = 0                           #[cm]
        self.zvel = zvel                        #[cm/ms]
        self.zpos_init = zpos_init              #[cm]

        # Now we check if any of those lists are of length > 1
        #if not self.any_to_fit():
            # We know all of these arrays will have EXACTLY 1 element
            # If they were empty, that would have been caught in arr_or_num()
            # If they had length > 1 then we would not be in this conditional
            # statement.
            # So we can set all of the current simulation values to the first
            # (and only) element of the corresponding arrays
        self.T1_f = self.T1_f_vals[0]
        self.T2_f = self.T2_f_vals[0]
        self.T1_s = self.T1_s_vals[0]
        self.alpha = self.alpha_vals[0]
        self.F = self.F_vals[0]
        self.ks = self.ks_vals[0]
        self.kf = self.kf_vals[0]
        self.CBV = self.CBV_vals[0]
        self.BAT = self.BAT_vals[0]
        # This way, we do not have to call iter() in any way to initialize
        # if the object doesnt have anything to iterate over anywqays.
        self.calc_R_T_vals()


    def __iter__(self):
        """
        This method makes the class Params iterable. It initializes the values and indices
        so that we may loop over the set of all parameters. It also sets two flags so that
        we know if we need to do any setup before simulating for each iteration.
        """

        # Initializing flags
        self.recompute_s = True     # Do not assume it has already been computed
        self.rescale_s = False      # We do not need to rescale if we just computed it

        # PARAMETER INDICES
        self.T1_f_ind = 0
        self.T2_f_ind = 0
        self.T1_s_ind = 0
        self.ks_ind = 0
        self.kf_ind = 0
        self.F_ind = 0
        self.CBV_ind = 0
        self.BAT_ind = 0
        self.alpha_ind = 0

        # INITIALIZING CURRENT VALUES
        self.T1_f = self.T1_f_vals[self.T1_f_ind]
        self.T2_f = self.T2_f_vals[self.T2_f_ind]
        self.T1_s = self.T1_s_vals[self.T1_s_ind]
        self.alpha = self.alpha_vals[self.alpha_ind]
        self.F = self.F_vals[self.F_ind]
        self.ks = self.ks_vals[self.ks_ind]
        self.kf = self.kf_vals[self.kf_ind]
        self.CBV = self.CBV_vals[self.CBV_ind]
        self.BAT = self.BAT_vals[self.BAT_ind]

        # Calculating apparent Rs and Ts
        self.calc_R_T_vals()

        return self


    def __next__(self):
        """
        This method moves along the iteration by one step. The function cascades
        through all fitting parameters to move forward as efficiently as possible.
        Once the first parameter on the list has done one full cycle, the next one
        is incremented by 1, and so on until we have made it all the way down the list.
        The parameters are incremented in the following order:

            --FIRST--
            - CBV
            - ks
            - kf
            - T1_f
            - T2_f
            - T1_s
            - F
            - alpha
            - BAT
            --LAST--

        The parameters that are incremented last involve the most setup, particulatly
        BAT, which requires that we recalculate a different s(t) signal from scratch.
        Parameters towards the beginning of the list require little to no setup for the
        MRFSim instance or the objects that it wraps.
        """
        # Update CBV index, if that index is != 0 then we do not update the next
        self.CBV_ind = (self.CBV_ind + 1) % np.size(self.CBV_vals)
        self.CBV = self.CBV_vals[self.CBV_ind]
        if self.CBV_ind: return None

        # Update ks index, if that index is != 0 then we do not update the next
        self.ks_ind = (self.ks_ind + 1) % np.size(self.ks_vals)
        self.ks = self.ks_vals[self.ks_ind]
        if self.ks_ind: return None

        # Update kf index, if that index is != 0 then we do not update the next
        self.kf_ind = (self.kf_ind + 1) % np.size(self.kf_vals)
        self.kf = self.kf_vals[self.kf_ind]
        if self.kf_ind: return None

        # Update T1_f index, if that index is != 0 then we do not update the next
        # We also recalculate the apparent R and T values
        self.T1_f_ind = (self.T1_f_ind + 1) % np.size(self.T1_f_vals)
        self.T1_f = self.T1_f_vals[self.T1_f_ind]
        if self.T1_f_ind: 
            self.calc_R_T_vals()
            return None

        # Update T2_f index, if that index is != 0 then we do not update the next
        # We also recalculate the apparent R and T values
        self.T2_f_ind = (self.T2_f_ind + 1) % np.size(self.T2_f_vals)
        self.T2_f = self.T2_f_vals[self.T2_f_ind]
        if self.T2_f_ind:
            self.calc_R_T_vals()
            return None
        
        # Update T1_s index, if that index is != 0 then we do not update the next
        # We also recalculate the apparent R and T values
        self.T1_s_ind = (self.T1_s_ind + 1) % np.size(self.T1_s_vals)
        self.T1_s = self.T1_s_vals[self.T1_s_ind]
        if self.T1_s_ind:
            self.calc_R_T_vals()
            return None
        
        # Update F index, if that index is != 0 then we do not update the next
        # We also recalculate the apparent R and T values
        self.F_ind = (self.F_ind + 1) % np.size(self.F_vals)
        self.F = self.F_vals[self.F_ind]
        if self.F_ind:
            self.calc_R_T_vals()
            self.rescale_s = True
            return None
        
        # Update alpha index, if that index is != 0 then we do not update the next
        # We also recalculate the apparent R and T values
        self.alpha_ind = (self.alpha_ind + 1) % np.size(self.alpha_vals)
        self.alpha = self.alpha_vals[self.alpha_ind]
        if self.alpha_ind:
            self.calc_R_T_vals()
            self.rescale_s = True
            return None
        
        # Update BAT index
        self.BAT_ind = (self.BAT_ind + 1) % np.size(self.BAT_vals)
        self.BAT = self.BAT_vals[self.BAT_ind]
        if self.BAT_ind == 0:
            # Once this happens, we have reached the end of the iteration
            raise StopIteration
            
        self.calc_R_T_vals()
        self.recompute_s = True
        return None
    

    def calc_R_T_vals(self):
        """
        This is a helper method that calculates apparent R and T values using 
        the current physiological values to be simulated. This method will
        be called after an update to T1_f, T2_f, T1_s, or F using the following
        formulas:

            R1_f_app = (F / lam) + (1 / T1_f)
            R2_f_app = (F / lam) + (1 / T2_f)

            T1_f_app = 1 / R1_f_app
            T2_f_app = 1 / R2_f_app

        These values are stored for convenience for the simulation to come.
        TODO: Eliminate R1s and T1s app (theyre not useful)
        """

        self.R1f_app = (self.F / self.lam) + (1 / self.T1_f)
        self.R2f_app = (self.F / self.lam) + (1 / self.T2_f)
        self.R1s_app = (1 / self.T1_s)

        self.T1f_app = 1 / self.R1f_app
        self.T2f_app = 1 / self.R2f_app


    def any_to_fit(self):
        """
        This is a helper method that checks if we have any parameters that
        will be fitted. In other words, this method checks if any of the 9 fitting
        dimensions have length > 1. If yes, then we will loop through at least one
        dimension of parameters and this function will return true, otherwise it will
        return false. It will also set a flag to reflect this so if this method is called
        twice we wont have to recompute.
        """

        # If we havent computed, we do it now
        if not hasattr(self, "fitting"):
            self.fitting = \
                (np.size(self.CBV_vals) > 1) | \
                (np.size(self.ks_vals) > 1) | \
                (np.size(self.kf_vals) > 1) | \
                (np.size(self.T1_f_vals) > 1) | \
                (np.size(self.T2_f_vals) > 1) | \
                (np.size(self.T1_s_vals) > 1) | \
                (np.size(self.F_vals) > 1) | \
                (np.size(self.alpha_vals) > 1) | \
                (np.size(self.BAT_vals) > 1)

        return self.fitting
    
    def print_inds(self):
        """
        Helper method that prints out the current indices of the parameters.
        """
        print("Current Indices: [ CBV:", self.CBV_ind, " , ks:", self.ks_ind, " , kf:",self.kf_ind, \
                " , T1_f:",self.T1_f_ind, " , T2_f:",self.T2_f_ind, " , T1_s:",self.T1_s_ind, " , F:", \
                self.F_ind, " , a:", self.alpha_ind, " , BAT:", self.BAT_ind, " ]")
        

            #     --FIRST--
            # - CBV
            # - ks
            # - kf
            # - T1_f
            # - T2_f
            # - T1_s
            # - F
            # - alpha
            # - BAT
            # --LAST--

    def get_shape(self):
        if not hasattr(self, "val_shape"):
            self.val_shape = ( \
            np.size(self.CBV_vals), \
            np.size(self.ks_vals), \
            np.size(self.kf_vals), \
            np.size(self.T1_f_vals), \
            np.size(self.T2_f_vals), \
            np.size(self.T1_s_vals), \
            np.size(self.F_vals), \
            np.size(self.alpha_vals), \
            np.size(self.BAT_vals), \
            )
        
        return self.val_shape
    

    def get_cur_idx(self):
        return ( \
            self.ks_ind, \
            self.CBV_ind, \
            self.kf_ind, \
            self.T1_f_ind, \
            self.T2_f_ind, \
            self.T1_s_ind, \
            self.F_ind, \
            self.alpha_ind, \
            self.BAT_ind, \
            )
    
    def get_num_combs(self):
        if not hasattr(self, "total_combs"):
            self.total_combs = np.prod(self.get_shape())

        return self.total_combs
    
    def get_comp_perc(self):
        return np.ravel_multi_index(self.get_cur_idx(), self.get_shape(), order="F") / self.get_num_combs()


        

def arr_or_num(arg):
    """
    Helper function that checks if the input is a numpy array, a number, or something else.
    This function will return a numpy array so that it can be iterated over in Params. If the
    input is already a numpy array, the function does nothing and returns it. If the input
    is a number, it wraps it in a numpy array (if this happens, we essentially arent iterating
    over that particular parameter). If the input is something else, the function raises a value
    error and stops execution.
    """
    if isinstance(arg, np.ndarray):
        # If the input is a numpy array, nothing to be done,
        # so we return it
        return arg
    elif isinstance(arg, Number):
        # If the input is a number, create an arr of length
        # 1 from it.
        return np.array([arg])
    else:
        raise ValueError("Error: Input arguments to the Params constructor must be numbers or ndarrays")