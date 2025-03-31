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
T = 1000                       # Simulation durration [ms]
dt = 10**-2                     # Simulation timestep [ms]
ntime = int(ceil(T / dt))       # Number of time samples
time = arange(ntime) * dt       # Time vector [ms]
time_s = time * 10**3           # Time vector for displaying in [s]

xpos = 0.0                      #[cm]
ypos = 0.0                      #[cm]
zvel = 10**-3                   #[cm/ms]
zpos_init = -T * zvel/2         #[cm]

d_psi = 0.05 * pi               # Phase gain of the RF pulse [rad]
psi_0 = 0.5 * pi                # Initial Phase of RF pulse [rad]
flip = 20                       # flip angle degrees

TR = 50                        # Repetition time for each pulse [ms]
PW = 10                         # width of each pulse [ms]

T1_b = 1600
T2_b = 150

rf_amp = 15* (10**-7)                # RF amplitude [T]
g_max = 0.00000006
g_ave = 0.000000006
max_slew = 0.5

# white matter values
F = 0.0
lam = 0.9
CBV = 0.005
BAT = 00
a_0 = 2

T1_s = 100


###########################



###########################
#   COORDINATE  INDICES   #
x = 0
y = 1
z = 2
n_dim = 3
###########################
