# Lorenz3D.py, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz3D
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

import sys
sys.path.append('../../')
import gbeespy as gbees # type: ignore

DIM_f = 3 # State dimension
DIM_h = 1 # Measurement dimension

# This function defines the dynamics model - required
def Lorenz3D(x, t, dx, coef):
    f1 = coef[0]*(x[1] - (x[0] + (dx[0]/2.0)))
    f2 = -(x[1] + (dx[1]/2.0)) - x[0]*x[2]
    f3 = -coef[1]*(x[2] + (dx[2]/2.0)) + x[0]*x[1] - coef[1]*coef[2]
    return [f1, f2, f3]

# This function defines the measurement model - required if MEASURE == True
def z(x, t, dx, coef):
    h1 = x[2]
    return [h1]

#==================================== Read in initial discrete measurement ==================================#
print("Reading in initial discrete measurement...\n")

P_DIR = "./results/python"      # Saved PDFs path
M_DIR = "./measurements"     # Measurement path
M_FILE = "measurement0.txt"  # Measurement file
M = gbees.Meas_create(DIM_f, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n")

dx = [None] * DIM_f                            # Grid width, default is half of the std. dev. from the initial measurement 
for i in range(DIM_f):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM_f, 5E-6, M.mean, dx) # Inputs: (dimension, probability threshold, center, grid width)    

coef = [4.0, 1.0, 48.0]                        # Lorenz3D trajectory attributes (sigma, beta, r)
T = gbees.Traj_create(len(coef), coef)         # Inputs: (# of coefficients, coefficients)

NUM_DIST = 5                                   # Number of distributions recorded per measurement
NUM_MEAS = 2                                   # Number of measurements
DEL_STEP = 20                                  # Number of steps per deletion procedure
OUTPUT_FREQ = 20                               # Number of steps per output to terminal
OUTPUT = True                                  # Write info to terminal
RECORD = True                                 # Write PDFs to .txt file
MEASURE = True                                 # Take discrete measurement updates
BOUNDS = False                                 # Add inadmissible regions to grid
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(Lorenz3D, z, None, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS)