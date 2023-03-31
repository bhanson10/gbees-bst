/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main() {

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    Grid G; G.T = 1; G.thresh = 0.00002; G.dt = 0.0005; G.start = {-11.5, -10, 9.5}; G.std = {1.0, 1.0, 1.0}; 
    G.del = {0.5,0.5,0.5}; G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2}; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; double per_cutoff = 42; int num_meas = 1; 
    int num_step = round(G.T/G.dt); int record_step = 1; int del_step = 17; 
    //====================================== End User Input =======================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,Lor); 
    
    std::cout << "Entering time marching..." << std::endl; 
    std::cout << std::endl; 

    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;
    double sim_time = 0; 
    std::string filename = "time_complexity_bst.txt"; std::ofstream myfile; myfile.open(filename);
    auto first_finish = start; 
    for(int nm = 0; nm < num_meas; nm++){
        D.Modify_pointset(G,Lor); 

        //Continous-time marching of PDF
        for (int i = 1; i <= num_step; i++){
            D.Modify_pointset(G,Lor); D.RHS_P(G, Lor); D.update_prob(G, D.P.root, D.a_count, D.tot_count); 
        
            if (i % record_step == 0){ 
                auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - first_finish; 
                myfile << D.tot_count << " " << elapsed.count() << std::endl;
                first_finish = finish; 
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
    myfile.close(); 

    return 0;
}