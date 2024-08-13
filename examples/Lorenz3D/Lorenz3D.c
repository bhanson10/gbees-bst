#include "../../gbees.c"

#define DIM 3

// This function defines the dynamics model - required
void Lorenz3D(double* x, double* dx, double* coef){
    double* v = (double*)malloc(DIM * sizeof(double)); 
    v[0] = coef[0]*(x[1]-(x[0]+(dx[0]/2.0)));
    v[1] = -(x[1]+(dx[1]/2.0))-x[0]*x[2];
    v[2] = -coef[1]*(x[2]+(dx[2]/2.0))+x[0]*x[1]-coef[1]*coef[2];
    for(int i = 0; i < DIM; i++){
        x[i] = v[i];
    }
    free(v); 
}

int main(){
    //=================================== Read in initial discrete measurement =================================//
    printf("\nReading in initial discrete measurement...\n\n");

    char* P_DIR = "<path_to_pdf>";     // Saved PDFs path
    char* M_DIR = ".";                 // Measurement path
    char* M_FILE = "/measurement.txt"; // Measurement file
    Meas M = Meas_create(DIM, M_DIR, M_FILE);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");

    double dx[DIM];                                           // Grid width, default is half of the std. dev. from the initial measurement 
    for(int i = 0; i < DIM; i ++){
        dx[i] = pow(M.cov[i][i],0.5)/2;
    }
    Grid G = Grid_create(DIM, 2E-5, M.mean, dx);              // Inputs: (dimension, probability threshold, center, grid width)       

    double coef[] = {4.0, 1.0, 48.0};                         // Lorenz3D trajectory attributes (sigma, beta, r)
    Traj T = Traj_create(3, coef); // Inputs: (# of coefficients, coefficients)

    int OUTPUT_FREQ = 20;                                     // Number of steps per output to terminal
    int DEL_STEP = 20;                                        // Number of steps per deletion procedure
    int NUM_DIST = 6;                                         // Number of distributions recorded per measurement
    int NUM_MEAS = 1;                                         // Number of measurements
    bool OUTPUT = true;                                       // Write info to terminal
    bool RECORD = true;                                       // Write PDFs to .txt file
    bool MEASURE = true;                                      // Take discrete measurement updates
    bool BOUNDS = false;                                       // Add inadmissible regions to grid
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(Lorenz3D, NULL, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE, BOUNDS);

    return 0;
}
