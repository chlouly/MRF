from globals import *


class Params:
    def __init__(self):
        # TODO: Get rid of values from the globals file
        T = 1000                       # Simulation durration [ms]
        dt = 10**-5                     # Simulation timestep [ms]
        ntime = int(ceil(T / dt))       # Number of time samples
        time = arange(ntime) * dt       # Time vector [ms]

        self.T1 = T1
        self.T2 = T2

        self.xpos = xpos                      #[cm]
        self.ypos = ypos                      #[cm]
        self.zvel = zvel                   #[cm/ms]
        self.zpos_init = zpos_init         #[cm]

        self.d_psi = d_psi               # Phase gain of the RF pulse [rad]
        self.psi_0 = psi_0                # Initial Phase of RF pulse [rad]

        self.TR = TR                        # Repetition time for each pulse [ms]
        self.PW = PW                        # width of each pulse [ms]

        self.rf_amp = rf_amp                # RF amplitude [T]
        self.g_max = g_max
        self.g_ave = g_ave
        self.max_slew = max_slew