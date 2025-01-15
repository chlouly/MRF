"""
All Global variables live in this file

Any simulation parameters should live here
"""

from numpy import ceil, arange

###########################
#        Constants        #
pi = 3.14159
gambar = 42570                  # Gyromagnetic coefficient [kHz/T]
gam = gambar * 2 * 3.14159      # Gamma [kRad/sT]
###########################



###########################
#  Experiment Parameters  #
T = 10000                       # Simulation durration [ms]
dt = 10**-2                     # Simulation timestep [ms]
ntime = int(ceil(T / dt))       # Number of time samples
time = arange(ntime) * dt       # Time vector [ms]
time_s = time * 10**3           # Time vector for displaying in [s]

xpos = 0.0                      #[cm]
ypos = 0.0                      #[cm]
zvel = 10**-3                   #[cm/ms]
zpos_init = -5.0                #[cm]

d_psi = 0.05 * pi               # Phase gain of the RF pulse [rad]
psi_0 = 0.5 * pi                # Initial Phase of RF pulse [rad]

TR = 500                        # Repetition time for each pulse [ms]
PW = 200                        # width of each pulse [ms]

rf_amp = 0.00005                # RF amplitude [T]
g_max = 0.0001
g_ave = 0.1
max_slew = 12000
rephase_fact = 4                # How many times longer the main puse is than the rephasing pulse
###########################



###########################
#   COORDINATE  INDICES   #
x = 0
y = 1
z = 2
n_dim = 3
###########################
