import sys
sys.path.append('../../')
import gbeespy as gbees

DIM = 3

# This function defines the dynamics model - required
def Lorenz3D(x, dx, coef):
    v1 = coef[0]*(x[1] - (x[0] + (dx[0]/2.0)))
    v2 = -(x[1] + (dx[1]/2.0)) - x[0]*x[2]
    v3 = -coef[1]*(x[2] + (dx[2]/2.0)) + x[0]*x[1] - coef[1]*coef[2]
    return [v1, v2, v3]

#==================================== Read in initial discrete measurement ==================================#
print("\nReading in initial discrete measurement...\n\n")

P_DIR = "./results/python"    # Saved PDFs path
M_DIR = "."                   # Measurement path
M_FILE = "/measurement0.txt"; # Measurement file
M = gbees.Meas_create(DIM, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n\n")

dx = [None] * DIM                             # Grid width, default is half of the std. dev. from the initial measurement 
for i in range(DIM):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM, 2E-5, M.mean, dx); # Inputs: (dimension, probability threshold, center, grid width)    

coef = [4.0, 1.0, 48.0]                       # Lorenz3D trajectory attributes (sigma, beta, r)
T = gbees.Traj_create(len(coef), coef);       # Inputs: (# of coefficients, coefficients)

OUTPUT_FREQ = 20;                             # Number of steps per output to terminal
DEL_STEP = 20;                                # Number of steps per deletion procedure
NUM_DIST = 6;                                 # Number of distributions recorded per measurement
NUM_MEAS = 1;                                 # Number of measurements
OUTPUT = True;                                # Write info to terminal
RECORD = True;                                # Write PDFs to .txt file
MEASURE = True;                               # Take discrete measurement updates
BOUNDS = False;                               # Add inadmissible regions to grid
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(Lorenz3D, None, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE, BOUNDS)