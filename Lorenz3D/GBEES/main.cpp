//
// Created by Benjamin Hanson on 5/7/24.
//

#include "BST.h"

int main(){

    //===================================== Begin User Input ====================================================
    std::cout << "Reading in user inputs..." << std::endl; std::cout << std::endl;

    const int NM = 1;                         // Number of measurements 
    std::string FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/Lorenz3D/Current/GBEES/cmake-build-debug/Data"; // Measurement file path
    std::string FILE_PATH_M = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/Lorenz3D/Current/GBEES/cmake-build-debug/Movie Data"; // Movie file path

    Grid G{};               // Grid object
    G.thresh = 2E-5;        // Probability threshold
    G.DIFF_B = false;       // Diffusion inclusion boolean
    G.diff = {0, 0, 0};     // Diffusion coefficient vector
    bool OUTPUT = true;     // Write info to terminal
    bool RECORD = true;     // Write PDFs to .txt file
    bool RECORD_F = false;  // Write frequent PDFs to .txt file (for movie)
    bool MEASURE = true;    // Take discrete measurement updates
    int output_freq = 20;   // Number of steps per output to terminal
    int record_f_freq = 20; // Number of steps per frequent record (for movie)
    int del_step = 25;      // Number of steps per deletion procedure
    int num_dist = 6;       // Number of distributions recorded per measurement
    //===========================================================================================================

    //===================================== Read in measurement/trajectory info =================================
    std::cout << "Reading in discrete measurements..." << std::endl; std::cout << std::endl;

    std::ifstream measurement_file(FILE_PATH + "/measurements.txt"); // Measurement file
    double means[NM][DIM];
    double stds[NM][DIM];
    int r = 0;
    std::string line;
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < NM)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            if(c < DIM){
                means[r][c] = std::stod(line);
            }else{
                stds[r][c-DIM] = std::stod(line);
            }
            c++;
        }
        r++;
    }

    G.epoch = {means[0][0], means[0][1], means[0][2]};  // Grid initial epoch
    G.del = {stds[0][0]/2, stds[0][1]/2, stds[0][2]/2}; // Grid width
    G.xh = {G.del[0]/2, G.del[1]/2, G.del[0]/3};        // Half grid width

    Traj Lor{}; // Trajectory object
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); Lor.sigma = std::stod(line); // Read in sigma
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); Lor.b = std::stod(line);  // Read in b
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); Lor.r = std::stod(line);  // Read in r
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); Lor.T = std::stod(line);  // Read in T

    double measure_time = Lor.T/NM;            // Time between measurements
    double record_time = measure_time/(num_dist-1); // Time between recording PDF

    Measurement m{};
    m.mean = {means[0][0], means[0][1], means[0][2]}; // Measurement mean
    m.std = {stds[0][0], stds[0][1], stds[0][2]};      // Measurement std

    measurement_file.close();
    //===========================================================================================================
    BST P;

    std::cout << "Initializing Distribution...\n" << std::endl;

    P.initialize_grid(G, Lor, m);     // Create initial GBEES distribution
    P.normalize_tree(G, P.root);   // Normalize distribution

    std::cout << "Entering time marching...\n" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed{};

    double tt = 0; int nm = 0; int record_count, step_count, record_f_count; double rt, mt; std::ofstream time_file;
    while(tt < Lor.T){ // Time of uncertainty propagation

        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; P.max_key(G, P.root, key); std::string filename;
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << tt;
        std::cout << " TU, Active/Total Cells: " << P.a_count << "/" << P.tot_count << ", Max key %: " << (double(key)/(pow(2,64)-1))*(100) << '%' << std::endl;
        if(RECORD){filename = FILE_PATH + "/M" + std::to_string(nm) + "/pdf_0.txt"; P.record_data(filename, G, 0);};
        if(RECORD_F){filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/pdf_0.txt"; P.record_data(filename, G, 0);};

        mt = 0; step_count = 0; record_count = 1; record_f_count = 1;
        while(mt < measure_time){ // Time between measurements 
            rt = 0;

            while(rt < record_time){ // Time between recording the PDF

                P.grow_tree_1(G,Lor);
                P.check_cfl_condition(G, P.root);
                G.dt = std::min(P.cfl_min_dt, record_time-rt); rt += G.dt;
                P.godunov_method(G);
                P.update_prob(G, P.root);
                P.normalize_tree(G, P.root);

                if (step_count % del_step == 0){ // Deletion procedure
                    P.prune_tree(G,P.root);
                }

                step_count+=1;
                if (step_count % output_freq == 0){ // Print size to terminal
                    finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
                    key = 0; P.max_key(G, P.root, key);
                    std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt;
                    std::cout << " TU, Active/Total Cells: " << P.a_count << "/" << P.tot_count << ", Max key %: " << (double(key)/(pow(2,64)-1))*(100) << '%' << std::endl;
                }

                if((RECORD_F)&&(step_count % record_f_freq == 0)){ // Record Movie data
                    filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/pdf_" + std::to_string(record_f_count) + ".txt"; P.record_data(filename, G, mt + rt);
                    record_f_count+=1;
                }

                P.cfl_min_dt = 1E10; time_file << std::setprecision(8) << mt + rt << std::endl;
            }
            if((OUTPUT)&&(step_count % output_freq != 0)){
                finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
                key = 0; P.max_key(G, P.root, key);
                std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt;
                std::cout << " TU, Active/Total Cells: " << P.a_count << "/" << P.tot_count << ", Max key %: " << (double(key)/(pow(2,64)-1))*(100) << '%' << std::endl;
            }

            if(RECORD){ // Record PDF 
                std::cout << std::endl;
                std::cout << "RECORDING PDF AT: " << tt + mt + rt << " TU..." << std::endl;
                std::cout << std::endl;
                filename = FILE_PATH + "/M" + std::to_string(nm) + "/pdf_" + std::to_string(record_count) + ".txt"; P.record_data(filename, G, mt + rt);
                record_count+=1;
            }
            if(RECORD_F){ // Record PDF 
                filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/pdf_" + std::to_string(record_count) + ".txt"; P.record_data(filename, G, mt + rt);
                record_f_count+=1;
            }

            mt += rt;
        }

        tt += mt;

        if((MEASURE)&&(tt < Lor.T)){ // Perform discrete measurement update
            std::cout << std::endl;
            std::cout << "PERFORMING BAYESIAN UPDATE AT: " << tt << " TU..." << std::endl;
            std::cout << std::endl;
            nm+=1;

            m.mean = {means[nm][0], means[nm][1], means[nm][2]};
            m.std = {stds[nm][0], stds[nm][1], stds[nm][2]};
            P.measurement_update(m, G, P.root); P.normalize_tree(G, P.root);
            P.prune_tree(G, P.root); // Deletion procedure
        }

        time_file.close(); measurement_file.close();
    }

    return 0;
}