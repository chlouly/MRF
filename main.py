from src import *
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
#from py_files = glob.glob("./src/*.py", recursive=True) + glob.glob("./src/*/*.py", recursive=True) 
print()
print()
print(os.getcwd())
print(py_files)
print()
print()
print()

extensions=[]
for filepath in py_files:
    # Convert file path to module name (e.g., src/myfastlib/core.pyx → myfastlib.core)
    module_name = os.path.splitext(filepath.replace("src/", "").replace(os.sep, "."))[0]
    extensions.append(Extension(name=module_name, sources=[filepath])).test.test_globals import *
py_files = glob.glob("./src/*.py", recursive=True) + glob.glob("./src/*/*.py", recursive=True) 
print()
print()
print(os.getcwd())
print(py_files)
print()
print()
print()

extensions=[]
for filepath in py_files:
    # Convert file path to module name (e.g., src/myfastlib/core.pyx → myfastlib.core)
    module_name = os.path.splitext(filepath.replace("src/", "").replace(os.sep, "."))[0]
    extensions.append(Extension(name=module_name, sources=[filepath]))
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

T1_f = 1000
T2_f = 150

rf_amp = 15* (10**-7)                # RF amplitude [T]
g_max = 0.00000006
g_ave = 0.000000006
max_slew = 0.5

# white py_files = glob.glob("./src/*.py", recursive=True) + glob.glob("./src/*/*.py", recursive=True) 
print()
print()
print(os.getcwd())
print(py_files)
print()
print()
print()

extensions=[]
for filepath in py_files:
    # Convert file path to module name (e.g., src/myfastlib/core.pyx → myfastlib.core)
    module_name = os.path.splitext(filepath.replace("src/", "").replace(os.sep, "."))[0]
    extensions.append(Extension(name=module_name, sources=[filepath]))matter values
F = 50 / 60 / 1000 / 100
lam = 0.9
CBV = 0.02
BAT = 4000
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

#import matplotlib.pyplot as plt

SCHED_PATH = "/home/clouly/DEV/sched_bldr"
DICT_PATH = "/home/clouly/DEV/dicts"


if __name__ == "__main__":
    # T1 and T2 fitting (generating dictionary)
    T1_f = np.linspace(300, 2000, 100)
    T2_f = np.linspace(100, 2000, 100)
    #T2_f = 1000
    flip = np.linspace(5, 25, 10)
    

    p = Params(T1_f, T2_f, T1_s, 0.0000, 0.0000, 0, lam, 0, 0, 0, BAT, 1, 1, flip)
    ps = MRFSim(p)

    ps.read_sched(SCHED_PATH + "/90006")

    ps.setup()

    ps.generate_dict(DICT_PATH + "/90006_GRE_BALL.h5")
