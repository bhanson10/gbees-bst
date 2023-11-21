/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main(){

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    std::ifstream measurement_file("./Data/measurements_SE.txt"); 
    double measurements[NM][DIM]; int r = 0; std::string line; std::istringstream iss;
    while(std::getline(measurement_file, line)) {
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')) {
            measurements[r][c] = std::stod(line); 
            c++;
        }
        r++;
    }
    std::ofstream time_file; 
    std::ofstream size_file; 
    //size_file.open("./MC/size_SE.txt");
    //===================================== Begin User Input ======================================
    Grid G; G.thresh = 1E-7; G.epoch = {measurements[0][0], measurements[0][1], measurements[0][2], measurements[0][3]};
    G.std = {1E-4, 1E-4, 5E-2, 5E-2}; G.del = {1E-4, 1E-4, 5E-2, 5E-2};
    G.pair = 2; /*1: Cantor, 2: Rosenberg*/ G.del_method = 1; /*1. Steps, 2. %*/
    Traj lyap; lyap.mu = 3.0542E-06; lyap.T = 0.03141584257;
    double measure_time = lyap.T; double record_time = measure_time/5; bool RECORD = true; bool MEASURE = true;
    int output_freq = 100; int del_step = 10;  
    //====================================== End User Input =======================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,lyap); D.normalize_tree(G, D.P.root);
    
    std::cout << "Entering time marching..." << std::endl; 
    std::cout << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;

    double tt = 0; int nm = 0; int record_count; double rt; double mt; int step_count; 
    while(tt < lyap.T){ // Time of uncertainty propagation

        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; D.max_key(G, D.P.root, key); 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << tt;
        std::cout << " TU, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key %: " << (key/(pow(2,64)-1))*(100) << '%' << std::endl; 
        std::string filename = "./Data/M" + std::to_string(nm) + "/pdf_0.txt"; D.Record_Data(filename, G, D.P.root, 0); 
        time_file.open("./Data/Times/time" + std::to_string(nm) + ".txt"); time_file << 0 << std::endl; 

        mt = 0; step_count = 0; record_count = 1;
        while(mt < measure_time){ // Time between measurements 
            rt = 0; 

            while(rt < record_time){ // Time between recording the PDF

                D.Modify_pointset(G,lyap); D.check_CFL_condition(G, D.P.root, D.cfl_min_dt); G.dt = std::min(D.cfl_min_dt, record_time-rt); rt += G.dt;  
                D.RHS(G, lyap); D.update_prob(G, D.P.root); D.normalize_tree(G, D.P.root);  
                
                step_count+=1; 
                if (step_count % output_freq == 0){ // Print size to terminal
                    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; uint64_t key = 0; D.max_key(G, D.P.root, key); 
                    std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt << " TU, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key %: " << (key/(pow(2,64)-1))*(100) << '%' << std::endl; 
                }    

                if (step_count % del_step == 0){ // Deletion procedure
                    D.P.root = D.prune_tree(G,D.P.root); D.Initialize_ik_nodes(G,D.P.root);
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
                std::string filename = "./Data/M" + std::to_string(nm) + "/pdf_" + std::to_string(record_count) + ".txt"; D.Record_Data(filename, G, D.P.root, mt + rt);
                record_count+=1;
            }

            mt += rt; 
        }

        tt += mt; 
        //size_file << D.a_count << std::endl;

        if((MEASURE)&&(tt < lyap.T)){ // Perform discrete measurement update
            std::cout << std::endl; 
            std::cout << "PERFORMING BAYESIAN UPDATE AT: " << tt << " TU..." << std::endl;
            std::cout << std::endl; 
            nm+=1; 

            Measurement m; 
            m.mean = {measurements[nm][0], measurements[nm][1], measurements[nm][2], measurements[nm][3]}; m.std = G.std; 
            D.measurement_update(m, G, D.P.root); D.normalize_tree(G, D.P.root);
            D.P.root = D.prune_tree(G, D.P.root); D.Initialize_ik_nodes(G, D.P.root);
        }
        
        time_file.close(); measurement_file.close(); size_file.close(); 
    }

    return 0;
}