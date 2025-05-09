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
zvel = 0                        #[cm/ms]
zpos_init = 0                   #[cm]

d_psi = 0.05 * pi               # Phase gain of the RF pulse [rad]
psi_0 = 0.5 * pi                # Initial Phase of RF pulse [rad]
flip = 20                       # flip angle degrees

TR = 50                        # Repetition time for each pulse [ms]
PW = 10                         # width of each pulse [ms]

T1_f = 1820
T2_f = 150

rf_amp = 15* (10**-7)                # RF amplitude [T]
g_max = 0.00000006
g_ave = 0.000000006
max_slew = 0.5

# white matter values
F = 0.00001
lam = 0.9
CBV = 0.005
BAT = 7000
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
