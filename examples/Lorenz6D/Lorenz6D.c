// Lorenz6D.c, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz6D
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include "../../gbees.h" 
#include "Lorenz6D.h"

// This function defines the dynamics model - required
void Lorenz6D(double* f, double* x, double t, double* dx, double* coef){
    f[0] = (x[1] - x[4]) * x[5] - x[0] + coef[0];
    f[1] = (x[2] - x[5]) * x[0] - x[1] + coef[0];
    f[2] = (x[3] - x[0]) * x[1] - x[2] + coef[0];
    f[3] = (x[4] - x[1]) * x[2] - x[3] + coef[0];
    f[4] = (x[5] - x[2]) * x[3] - x[4] + coef[0];
    f[5] = (x[0] - x[3]) * x[4] - x[5] + coef[0];
}

int main(void){
    //=================================== Read in initial discrete measurement =================================//
    printf("Reading in initial discrete measurement...\n\n");

    char* P_DIR = "./results/c";     // Saved PDFs path
    char* M_DIR = "./measurements";    // Measurement path
    char* M_FILE = "measurement0.txt"; // Measurement file
    Meas M = Meas_create(DIM_f, M_DIR, M_FILE);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");

    double dx[DIM_f];                              // Grid width, default is half of the std. dev. from the initial measurement 
    for(int i = 0; i < DIM_f; i ++){
        dx[i] = pow(M.cov[i][i],0.5)/2;
    }
    Grid G = Grid_create(DIM_f, 8E-9, M.mean, dx); // Inputs: (dimension, probability threshold, center, grid width)       

    double coef[] = {4.0};              // Lorenz6D trajectory attributes (sigma, beta, r)
    Traj T = Traj_create(1, coef);                 // Inputs: (# of coefficients, coefficients)

    int NUM_DIST = 2;                              // Number of distributions recorded per measurement
    int NUM_MEAS = 1;                              // Number of measurements
    int DEL_STEP = 20;                             // Number of steps per deletion procedure
    int OUTPUT_FREQ = 20;                          // Number of steps per output to terminal
    bool OUTPUT = true;                            // Write info to terminal
    bool RECORD = true;                           // Write PDFs to .txt file
    bool MEASURE = false;                           // Take discrete measurement updates
    bool BOUNDS = false;                           // Add inadmissible regions to grid
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(Lorenz6D, NULL, NULL, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM_f, OUTPUT, RECORD, MEASURE, BOUNDS);

    return 0;
}
