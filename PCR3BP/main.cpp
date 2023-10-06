/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main(){

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    double LU = 389703; double TU = 382981; 
    Grid G; /*G.T = 6.5293315073396165;*/ G.T = 3; G.thresh = 1E-5; 
    G.epoch = {.64705017976217949, 0, 3.0487778939009659e-13, .75905615964008910};
    G.std = {5/LU, 5/LU, (5e-5)*(TU/LU),(5e-5)*(TU/LU)}; G.del = {2E-3, 2E-3, 4E-3, 4E-3}; 
    G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2, G.del[3]/2}; G.rk = 1; /*1: EE, 2: RK2, 3: RK3*/
    G.pair = 2; /*1: Cantor, 2: Rosenberg*/ G.del_method = 1; /*1. Steps, 2. %*/
    Traj lyap; lyap.mu = 1.215058560962404e-2; double per_cutoff = 42; int num_meas = 1;
    int output_freq = 100; int record_freq = 1; int del_step = 20; 
    //====================================== End User Input =======================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,lyap); D.Modify_pointset(G,lyap); D.get_size(G, D.P.root, D.a_count, D.tot_count, D.max_key); 
    
    std::cout << "Entering time marching..." << std::endl; 
    std::cout << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;
    double sim_time = 0; 

    for(int nm = 0; nm < 1; nm++){
        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << sim_time;
        std::cout << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key overflow: " << (D.max_key/(pow(2,64)-1))*(100) << '%' << std::endl; 
        std::string filename = "./Data/pdf_" + std::to_string(nm) + "-0.txt"; D.Record_Data(filename, G, D.P.root, 0); 
        D.a_count = 0; D.tot_count = 1; D.max_key = 0;

        //Continous-time marching of PDF
        int step_count = 1; int record_count = 0; 
        while (sim_time < G.T){            
            D.Modify_pointset(G,lyap); D.check_CFL_condition(G, D.P.root, D.cfl_min_dt); G.dt = std::min(D.cfl_min_dt, G.T-sim_time); sim_time += G.dt;  
            
            // RK Update
            for(int rk = 1; rk <= G.rk; rk++){
                D.RHS_P(G, lyap, rk); D.update_prob(G, D.P.root, rk); 
            }
        
             // Recording PDF
            if (step_count % record_freq == 0){ 
                std::string filename = "./Data/pdf_" + std::to_string(nm) + "-" + std::to_string(record_count) + ".txt"; D.Record_Data(filename, G, D.P.root, sim_time);
                record_count+=1; 
            }

            // Output Size
            if (step_count % output_freq == 0){ 
                auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;  D.get_size(G, D.P.root, D.a_count, D.tot_count, D.max_key); 
                std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << sim_time << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key overflow: " << (D.max_key/(pow(2,64)-1))*(100) << '%' << std::endl; 
                D.a_count = 0; D.tot_count = 1; 
            }    

            // Deletion Steps
            if (step_count % del_step == 0){
                D.P.root = D.delete_neighbors(G,D.P.root); 
                D.Initialize_ik_nodes(G,D.P.root);
            } 
            
            /*
            int i_count = D.tot_count - D.a_count; double i_per = ((double)i_count / (double)D.tot_count)*100; 
            if (i_per >= per_cutoff){
                D.P.root = D.delete_neighbors(G,D.P.root); 
                D.Initialize_ik_nodes(G,D.P.root); 
            }
            */

            D.cfl_min_dt = 1E10; 
            step_count+=1; 
        }   
        finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;  D.get_size(G, D.P.root, D.a_count, D.tot_count, D.max_key); 
        std::cout << "Timestep: " << nm << "-" << step_count-1 << ", Program time: " << elapsed.count() << " s, Sim. time: " << sim_time << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << std::endl;

        /*Discrete-time update of PDF
        std::cout << std::endl; 
        std::cout << "Receieved discrete-time measurement. Performing Bayesian update..." << std::endl;
        std::cout << std::endl; 

        Measurement m; m.mean[2] = -15; m.unc[2] = 2; 
        double C = 0; D.measurement_update(m, G, D.P.root, C); D.normalize_prob(D.P.root, C);
        D.P.root = D.delete_neighbors(G,D.P.root); D.Initialize_ik_nodes(G,D.P.root); 
        */
    }
    return 0;
}