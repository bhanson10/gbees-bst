/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main(){

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    double LU = 389703; double TU = 382981; 
    Grid G; /*G.T = 2.1728022077038327E;*/ G.T = 1E-4; G.thresh = 1E-8; G.dt = 1E-5; 
    G.epoch = {8.7635429542029941E-1, 5.1148290147058553E-27, -1.9192404368074575E-1, -4.9903376573332209E-14, 2.3007020677463236E-1, -1.4808938468336447E-13};
    G.std = {5E-3, 5E-3, 5E-3, 1E-8, 1E-8, 1E-8}; G.del = G.std; G.rk = 1; /*1: EE, 2: RK2, 3: RK3*/
    G.pair = 2; /*1: Cantor, 2: Rosenberg*/ G.del_method = 1; /*1. Steps, 2. %*/
    Traj lyap; lyap.mu = 1.215058560962404e-2; double per_cutoff = 42; int num_meas = 1; int del_step = 20; 
    int num_step = round(G.T/G.dt); int record_step = round(G.T/(100*G.dt));
    //====================================== End User Input =======================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,lyap); 
    
    std::cout << "Entering time marching..." << std::endl; 
    std::cout << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;
    double sim_time = 0; 

    for(int nm = 0; nm < 1; nm++){
        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
        std::string filename = "./Data/pdf_" + std::to_string(nm) + "-0.txt";
        D.Modify_pointset(G,lyap); D.Record_Data(filename, G, D.P.root); D.get_size(G, D.P.root, D.a_count, D.tot_count, D.max_key); 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << sim_time;
        std::cout << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key overflow: " << (D.max_key/(pow(2,64)-1))*(100) << '%' << std::endl; 
        D.a_count = 0; D.tot_count = 1; D.max_key = 0; 
    
        //Continous-time marching of PDF
        for (int i = 1; i <= num_step; i++){
            D.Modify_pointset(G,lyap); 
            
            // RK Update
            for(int rk = 1; rk <= G.rk; rk++){
                D.RHS_P(G, lyap, rk); D.update_prob(G, D.P.root, rk); 
            }
        
            // Recording PDF
            if (i % record_step == 0){ 
                auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;  D.get_size(G, D.P.root, D.a_count, D.tot_count, D.max_key); 
                std::cout << "Timestep: " << nm << "-" << i << ", Program time: " << elapsed.count() << " s, Sim. time: " << sim_time + i*G.dt;
                std::cout  << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << ", Max key overflow: " << (D.max_key/(pow(2,64)-1))*(100) << '%' << std::endl; 
                std::string filename = "./Data/pdf_" + std::to_string(nm) + "-" + std::to_string(i) + ".txt"; D.Record_Data(filename, G, D.P.root);
            }    
            // Deletion Steps
            if(G.del_method==1){
                if (i % del_step == 0){
                    D.P.root = D.delete_neighbors(G,D.P.root); 
                    D.Initialize_ik_nodes(G,D.P.root);
                } 
            }else if(G.del_method==2){
                int i_count = D.tot_count - D.a_count; double i_per = ((double)i_count / (double)D.tot_count)*100; 
                if (i_per >= per_cutoff){
                    D.P.root = D.delete_neighbors(G,D.P.root); 
                    D.Initialize_ik_nodes(G,D.P.root); 
                }
            }
            D.a_count = 0; D.tot_count = 1; D.max_key = 0; 
        }   
        sim_time += num_step*G.dt;  

        //Discrete-time update of PDF
        std::cout << std::endl; 
        std::cout << "Receieved discrete-time measurement. Performing Bayesian update..." << std::endl;
        std::cout << std::endl; 

        Measurement m; m.mean[2] = -15; m.unc[2] = 2; 
        double C = 0; D.measurement_update(m, G, D.P.root, C); D.normalize_prob(D.P.root, C);
        D.P.root = D.delete_neighbors(G,D.P.root); D.Initialize_ik_nodes(G,D.P.root); 
    }
    return 0;
}