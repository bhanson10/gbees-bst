/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main(){

    //===================================== Begin User Input ====================================================
    std::cout << "Reading in user inputs..." << std::endl; std::cout << std::endl; 

    const int NM = 1;                                   // Number of measurements 
    std::string FILE_PATH = "./Data/Sun-Earth";         // Measurement file path
    std::string FILE_PATH_M = "./Movie Data/Sun-Earth"; // Movie file path

    Grid G;                 // Grid object
    G.thresh = 1.5E-7;        // Probability threshold
    G.pair = 2;             // 1: Cantor, 2: Rosenberg
    G.DIFF_B = true;        // Diffusion inclusion boolean
    G.diff = {1E-5, 1E-5, 1E-5, 1E-5};  // Diffusion coefficient vector
    bool RECORD = true;     // Write PDFs to .txt file
    bool RECORD_F = false;   // Write frequent PDFs to .txt file (for movie)
    bool MEASURE = true;    // Take discrete measurement updates
    int output_freq = 100;  // Number of steps per output to terminal
    int record_f_freq = 10; // Number of steps per frequent record
    int del_step = 10;      // Number of steps per deletion procedure
    int num_dist = 6;       // Number of distributions recorded per measurement
    //===========================================================================================================

    //===================================== Read in measurement/trajectory info =================================
    std::cout << "Reading in discrete measurements..." << std::endl; std::cout << std::endl; 

    std::ifstream measurement_file(FILE_PATH + "/measurements.txt"); // Measurement file
    double means[NM][DIM]; 
    double stds[NM][DIM];
    int r = 0; 
    std::string line; 
    std::istringstream iss; 
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

    G.epoch = {means[0][0], means[0][1], means[0][2], means[0][3]}; // Grid initial epoch
    G.del = {stds[0][0], stds[0][1], stds[0][2], stds[0][3]};       // Grid width
  
    Traj lyap; // Trajectory object
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.mu = std::stod(line); // Read in mu
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.T = std::stod(line);  // Read in T

    double measure_time = lyap.T/NM;            // Time between measurements
    double record_time = measure_time/num_dist; // Time between recording PDF

    Measurement m; 
    m.mean = {means[0][0], means[0][1], means[0][2], means[0][3]}; // Measurement mean
    m.std = {stds[0][0], stds[0][1], stds[0][2], stds[0][3]};      // Measurement std

    measurement_file.close();
    //===========================================================================================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; std::cout << std::endl; 
    
    D.Initialize_D(G,lyap, m);     // Create initial GBEES distribtion
    D.normalize_tree(G, D.P.root); // Normalize distribution
    
    std::cout << "Entering time marching..." << std::endl; std::cout << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); // Timing object
    std::chrono::duration<double> elapsed;

    double tt = 0; int nm = 0; int record_count, step_count, record_f_count; double rt, mt; std::ofstream time_file; 
    while(tt < lyap.T){ // Time of uncertainty propagation

        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; D.max_key(G, D.P.root, key); 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << tt;
        std::cout << " TU, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key %: " << (key/(pow(2,64)-1))*(100) << '%' << std::endl; 
        std::string filename = FILE_PATH + "/M" + std::to_string(nm) + "/pdf_0.txt"; D.Record_Data(filename, G, D.P.root, 0); 
        if(RECORD_F){filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/pdf_0.txt"; D.Record_Data(filename, G, D.P.root, 0);};
        time_file.open(FILE_PATH + "/Times/time" + std::to_string(nm) + ".txt"); time_file << 0 << std::endl; 

        mt = 0; step_count = 0; record_count = 1; record_f_count = 1; 
        while(mt < measure_time){ // Time between measurements 
            rt = 0; 

            while(rt < record_time){ // Time between recording the PDF

                D.Modify_pointset(G,lyap); D.check_CFL_condition(G, D.P.root, D.cfl_min_dt); G.dt = std::min(D.cfl_min_dt, record_time-rt); rt += G.dt;  
                D.RHS(G, lyap); D.update_prob(G, D.P.root); D.normalize_tree(G, D.P.root);  

                if (step_count % del_step == 0){ // Deletion procedure
                    D.P.root = D.prune_tree(G, lyap, D.P.root);
                } 
                
                step_count+=1; 
                if (step_count % output_freq == 0){ // Print size to terminal
                    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; D.max_key(G, D.P.root, key); 
                    std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt << " TU, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key %: " << (key/(pow(2,64)-1))*(100) << '%' << std::endl; 
                }    

                if(RECORD_F){ // Record Movie data
                    if (step_count % record_f_freq == 0){ 
                        std::string filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/pdf_" + std::to_string(record_f_count) + ".txt"; D.Record_Data(filename, G, D.P.root, mt + rt);
                        record_f_count+=1;
                    }
                }

                D.cfl_min_dt = 1E10; time_file << std::setprecision(8) << mt + rt << std::endl; 
            }
            if(step_count % output_freq != 0){
                finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; D.max_key(G, D.P.root, key); 
                std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt << " TU, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key %: " << (key/(pow(2,64)-1))*(100) << '%' << std::endl; 
            }
 
            if(RECORD){ // Record PDF 
                std::cout << std::endl; 
                std::cout << "RECORDING PDF AT: " << tt + mt + rt << " TU..." << std::endl;
                std::cout << std::endl; 
                std::string filename = FILE_PATH + "/M" + std::to_string(nm) + "/pdf_" + std::to_string(record_count) + ".txt"; D.Record_Data(filename, G, D.P.root, mt + rt);
                record_count+=1;
            }

            mt += rt; 
        }

        tt += mt; 

        if((MEASURE)&&(tt < lyap.T)){ // Perform discrete measurement update
            std::cout << std::endl; 
            std::cout << "PERFORMING BAYESIAN UPDATE AT: " << tt << " TU..." << std::endl;
            std::cout << std::endl; 
            nm+=1; 

            Measurement m; 
            m.mean = {means[nm][0], means[nm][1], means[nm][2], means[nm][3]}; 
            m.std = {stds[nm][0], stds[nm][1], stds[nm][2], stds[nm][3]}; 
            D.measurement_update(m, G, D.P.root); D.normalize_tree(G, D.P.root);
            D.P.root = D.prune_tree(G, lyap, D.P.root); // Deletion procedure
        }
        
        time_file.close(); 
    }
    return 0;
}