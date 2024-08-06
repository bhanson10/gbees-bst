import gbeespy as gbees

DIM = 6

def CR3BP(x, dx, T):
    r1 = ((x[0] + T.coef[0])**2 + (x[1])**2 + (x[2])**2)**1.5
    r2 = ((x[0] - 1 + T.coef[0])**2 + (x[1])**2 + (x[2])**2)**1.5
    v1 = x[3]
    v2 = x[4]
    v3 = x[5]
    v4 = 2*x[3]+x[0]-(T.coef[0]*(x[0]-1+T.coef[0])/r2)-((1-T.coef[0])*(x[0]+T.coef[0])/r1)
    v5 = -2*x[2]+x[1]-(T.coef[0]*x[1]/r2)-((1-T.coef[0])*x[1]/r1)
    v6 = -(T.coef[0]*x[2]/r2)-((1-T.coef[0])*x[2]/r1)
    return [v1, v2, v3, v4, v5, v6]

#==================================== Read in initial discrete measurement ==================================#
print("\nReading in initial discrete measurement...\n\n")

P_DIR = "./Data/PCR3BP/PDFs"
M_DIR = "./Data/PCR3BP/Measurements";
M_FILE = "/measurement0.txt"; 
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

OUTPUT = True;                                # Write info to terminal
RECORD = True;                                # Write PDFs to .txt file
MEASURE = True;                               # Take discrete measurement updates
OUTPUT_FREQ = 20;                             # Number of steps per output to terminal
DEL_STEP = 20;                                # Number of steps per deletion procedure
NUM_DIST = 17;                                # Number of distributions recorded per measurement
NUM_MEAS = 1;                                 # Number of measurements
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(CR3BP, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE)