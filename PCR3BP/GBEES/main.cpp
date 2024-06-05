//
// Created by Benjamin Hanson on 5/7/24.
//

#include "BST.h"

int main(){

    //===================================== Begin User Input ====================================================
    std::cout << "Reading in user inputs..." << std::endl; std::cout << std::endl;

    Grid G(DIM);        // Grid object
    G.thresh = 1E-7;       // Probability threshold
    G.DIFF_B = false;      // Diffusion inclusion boolean
    G.diff = {0, 0, 0, 0}; // Diffusion coefficient vector
    bool OUTPUT = true;    // Write info to terminal
    bool RECORD = true;    // Write PDFs to .txt file
    int del_step = 10;     // Number of steps per deletion procedure
    int num_dist = 16;     // Number of distributions recorded per revolution
    int body = 0;          // 0: Europa
    int mn = 0;            // #: M#

    std::string FILE_PATH;
    if (body == 0) {
        FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/PCR3BP/GBEES/cmake-build-debug/Epochs/Europa";
    }
    //===========================================================================================================

    //===================================== Read in measurement/trajectory info =================================
    std::cout << "Reading in discrete measurements..." << std::endl; std::cout << std::endl;

    std::ifstream measurement_file(FILE_PATH + "/measurements.txt"); // Measurement file
    double means[1][DIM];
    int r = 0;
    std::string line;
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < 1)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            means[r][c] = std::stod(line);
            c++;
        }
        r++;
    }

    r = 0;
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < DIM)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            G.cov(r,c) = std::stod(line);
            c++;
        }
        r++;
    }

    G.epoch = {means[0][0], means[0][1], means[0][2], means[0][3]}; // Grid initial epoch
    G.del   = {pow(G.cov(0,0),0.5)/2, pow(G.cov(1,1),0.5)/2, pow(G.cov(2,2),0.5)/2, pow(G.cov(3,3),0.5)/2}; // Grid width

    Traj lyap{}; // Trajectory object
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.mu = std::stod(line); // Read in mu
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.T = std::stod(line);  // Read in T

    double record_time = lyap.T/(num_dist); // Time between recording PDF

    Measurement m{};
    m.mean = {means[0][0], means[0][1], means[0][2], means[0][3]}; // Measurement mean
    m.std  = {pow(G.cov(0,0),0.5), pow(G.cov(1,1),0.5), pow(G.cov(2,2),0.5), pow(G.cov(3,3),0.5)}; // Measurement std

    measurement_file.close();
    //===========================================================================================================
    BST P;

    std::cout << "Initializing Distribution...\n" << std::endl;

    P.initialize_grid(G, lyap, m);  // Create initial GBEES distribution
    P.normalize_tree(G, P.root); // Normalize distribution
    P.prune_tree(G,P.root);      // Prune tree

    std::cout << "Entering time marching...\n" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed{};

    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; P.max_key(G, P.root, key);
    std::cout << "Program time: " << elapsed.count() << " s, Sim. time: " << 0 << " TU";
    std::cout << ", Active/Total Cells: " << P.a_count << "/" << P.tot_count << ", Max key %: " << (double(key)/(pow(2,64)-1))*(100) << '%' << std::endl;
    if(RECORD){std::string filename = FILE_PATH + "/M" + std::to_string(mn) + "/pdf_0.txt";  P.record_data(filename, G, 0);};

    double tt = 0; double rt; int record_count = 1; int step_count = 0;
    while(tt < lyap.T){ // Time of uncertainty propagation
        rt = 0;
        while(rt < record_time){ // Time between recording the PDF

            P.grow_tree_1(G,lyap);
            P.check_cfl_condition(G, P.root);
            G.dt = std::min(P.cfl_min_dt, record_time-rt); rt += G.dt;
            P.godunov_method(G);
            P.update_prob(G, P.root);
            P.normalize_tree(G, P.root);

            if (step_count % del_step == 0){ // Deletion procedure
                P.prune_tree(G,P.root);
            }

            P.cfl_min_dt = 1E10;
        }

        tt += rt;

        if(OUTPUT){
            finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
            key = 0; P.max_key(G, P.root, key);
            std::cout << "Program time: " << elapsed.count() << " s, Sim. time: " << tt << " TU";
            std::cout << ", Active/Total Cells: " << P.a_count << "/" << P.tot_count << ", Max key %: " << (double(key)/(pow(2,64)-1))*(100) << '%' << std::endl;
        }

        if(RECORD){ // Record PF Scatter Plot
            std::string filename = FILE_PATH + "/M" + std::to_string(mn) + "/pdf_" + std::to_string(record_count) + ".txt"; P.record_data(filename, G, tt);
            record_count+=1;
        }
    }

    return 0;
}