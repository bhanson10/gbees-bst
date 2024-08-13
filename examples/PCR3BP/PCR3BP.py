import sys
sys.path.append('../../')
import gbeespy as gbees  # type: ignore

DIM = 4

# This function defines the dynamics model - required
def PCR3BP(x, dx, coef):
    r1 = ((x[0] + coef[0])**2 + (x[1])**2)**1.5
    r2 = ((x[0] - 1 + coef[0])**2 + (x[1])**2)**1.5
    v1 = x[2]
    v2 = x[3]
    v3 = 2*x[3]+x[0]-(coef[0]*(x[0]-1+coef[0])/r2)-((1-coef[0])*(x[0]+coef[0])/r1)
    v4 = -2*x[2]+x[1]-(coef[0]*x[1]/r2)-((1-coef[0])*x[1]/r1)
    return [v1, v2, v3, v4]

# This function defines the initial grid boundaries - optional
def PCR3BP_J(x, coef):
    r1 = ((x[0]+coef[0])**2+(x[1])**2)**0.5
    r2 = ((x[0]-1+coef[0])**2+(x[1])**2)**0.5
    J = (x[0])**2.0 + (x[1])**2.0 + (2*(1-coef[0])/r1) + (2*coef[0]/r2) + coef[0]*(1 - coef[0]) - ((x[2])**2.0 + (x[3])**2.0)
    return J

#==================================== Read in initial discrete measurement ==================================#
print("\nReading in initial discrete measurement...\n\n")

P_DIR = "<path_to_pdf>"      # Saved PDFs path
M_DIR = "."                  # Measurement path
M_FILE = "/measurement.txt"; # Measurement file
M = gbees.Meas_create(DIM, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n\n")

dx = [None] * DIM                             # Grid width, default is half of the std. dev. from the initial measurement 
for i in range(DIM):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM, 1E-7, M.mean, dx); # Inputs: (dimension, probability threshold, center, grid width)    
 
coef = [2.528017528540000E-5]                 # PCR3BP trajectory attributes (mu)
T = gbees.Traj_create(len(coef), coef);       # Inputs: (# of coefficients, coefficients)

OUTPUT_FREQ = 20;                             # Number of steps per output to terminal
DEL_STEP = 20;                                # Number of steps per deletion procedure
NUM_DIST = 17;                                # Number of distributions recorded per measurement
NUM_MEAS = 1;                                 # Number of measurements
OUTPUT = False;                               # Write info to terminal
RECORD = True;                                # Write PDFs to .txt file
MEASURE = True;                               # Take discrete measurement updates
BOUNDS = True;                                # Add inadmissible regions to grid
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(PCR3BP, PCR3BP_J, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE, BOUNDS)