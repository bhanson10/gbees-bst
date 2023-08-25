/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main() {

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    Grid G; G.T = 1; G.thresh = 0.00002; G.dt = 0.0005; G.epoch = {-11.5, -10, 9.5}; G.std = {1.0, 1.0, 1.0}; 
    G.del = {0.5,0.5,0.5}; G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2}; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; double per_cutoff = 42; int num_meas = 1; 
    int num_step = round(G.T/G.dt); int record_step = round(G.T/(5*G.dt)); int del_step = 26; 
    //====================================== End User Input =======================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,Lor); 
    
    std::cout << "Entering time marching..." << std::endl; 
    std::cout << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;
    double sim_time = 0; 
    
    for(int nm = 0; nm < num_meas; nm++){
        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
        std::string filename = "pdf_" + std::to_string(nm) + "-0.txt";
        D.Modify_pointset(G,Lor); D.Record_Data(filename, G, D.P.root); D.get_size(G, D.P.root, D.a_count, D.tot_count); 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << sim_time << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << std::endl;
        D.a_count = 0; D.tot_count = 1; 

        //Continous-time marching of PDF
        for (int i = 1; i <= num_step; i++){
            D.Modify_pointset(G,Lor); 
            
            for(int rk = 0; rk < 2; rk++){
                D.RHS_P(G, Lor, rk); D.update_prob(G, D.P.root, rk); 
            }
        
            if (i % record_step == 0){ 
                auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;  D.get_size(G, D.P.root, D.a_count, D.tot_count); 
                std::cout << "Timestep: " << nm << "-" << i << ", Program time: " << elapsed.count() << " s, Sim. time: " << sim_time + i*G.dt << " s, Active/Total Cells: " << D.a_count << "/" << D.tot_count << std::endl;
                std::string filename = "pdf_" + std::to_string(nm) + "-" + std::to_string(i) + ".txt"; D.Record_Data(filename, G, D.P.root);
            }    

            /*
            int i_count = D.tot_count - D.a_count; double i_per = ((double)i_count / (double)D.tot_count)*100; 
            if (i_per >= per_cutoff){
                D.P.root = D.delete_neighbors(G,D.P.root); 
                D.Initialize_ik_nodes(G,D.P.root); 
            }
            */
            
            if (i % del_step == 0){
                D.P.root = D.delete_neighbors(G,D.P.root); 
                D.Initialize_ik_nodes(G,D.P.root); 
            }

            D.a_count = 0; D.tot_count = 1; 
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