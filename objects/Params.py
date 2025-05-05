##########################################################################
#   This file contains the class definition for the Params object        #
#                                                                        #
#   Code written by Christopher Louly (clouly@umich.edu) 2025            #
##########################################################################


class Params:
    """
    This class holds all simulation parameters that are not related layout and timing
    of the pulse sequence. 
    """


    def __init__(self, T1_f, T2_f, T1_s, ks, kf, F, lam, zvel, zpos_init, CBV, BAT, M0_f, M0_s, alpha = 0.86):
        """
        This method initializes an instance of the Params object.

        Input Parameters:
            T1_f:       Free water (tissue) T1 [ms]
            T2_f:       Free water (tissue) T2 [ms]
            T1_s:       Semisolid pool (spins bound by macro-molecules) T1 [ms]
            ks:         Magnetization transfer rate (Semisolid -> Free Water) [ms]
            kf:         Magnetization transfer rate (Free Water -> Semisolid) [ms]
            F:          Flow Rate [mL / g / ms] ~= [1 / ms]
            lam:        Blood Brain Barrier Coeff [mL / g]
            zvel:       Observation Spin z velocity [cm / ms]
            zpos_init:  Observation spin initial z position [cm]
            CBV:        Cerebral Blood Volume [unitless]
            BAT:        Bollus Arrival Time [ms]
            M0_f:       Longitudinal Equilibrium Magnetization of the free water pool
            M0_s:       Longitudinal Equilibrium Magnetization of the semisolid pool
            alpha:      Labelling efficiency [unitless]
        """
        self.T1_f = T1_f
        self.T2_f = T2_f
        self.T1_s = T1_s
        self.T1_b = 1600        # Blood T1 (set to the typical 1600 ms)   

        self.alpha = alpha

        self.F = F
        self.lam = lam

        self.R1f_app = (self.F / self.lam) + (1 / self.T1_f)
        self.R2f_app = (self.F / self.lam) + (1 / self.T2_f)
        self.R1s_app = (1 / self.T1_s)

        self.T1f_app = 1 / self.R1f_app
        self.T2f_app = 1 / self.R2f_app

        self.ks = ks
        self.kf = kf

        self.M0_f = M0_f
        self.M0_s = M0_s

        self.f = M0_s / M0_f

        self.xpos = 0                           #[cm]
        self.ypos = 0                           #[cm]
        self.zvel = zvel                        #[cm/ms]
        self.zpos_init = zpos_init              #[cm]

        self.CBV = CBV
        self.BAT = BAT