/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main() {

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    Grid G; G.T = 1; G.thresh = 0.00002; G.dt = 0.0005; G.start = {-11.5, -10, 9.5}; G.std = {1.0, 1.0, 1.0}; 
    G.del = {0.5,0.5,0.5}; G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2}; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; 
    int num_step = round(G.T/G.dt); 
    //====================================== End User Input =======================================
    const int MAX_PER = 60; 
    const int START_PER = 30; 
    const int DEL_PER = 2; 
    const int ITERATIONS = 3; 
    const int LENGTH = ((MAX_PER-START_PER)/DEL_PER)+1;
    int count = 0; std::array<double, LENGTH> time; 
    
     for (int a = START_PER; a <= MAX_PER; a+=DEL_PER){
        double avg_time = 0; 
        std::cout << "Percentage/deletion: " << a << "%" << std::endl; 
        for (int b = 1; b <= ITERATIONS; b++){
            GBEES D;
            
            D.Initialize_D(G,Lor); D.Modify_pointset(G,Lor); 
            
            auto start = std::chrono::high_resolution_clock::now(); 
            std::chrono::duration<double> elapsed;

            for (int i = 1; i <= num_step; i++){
                D.Modify_pointset(G,Lor); D.RHS_P(G, Lor); D.update_prob(G, D.P.root, D.a_count, D.tot_count); 

                int i_count = D.tot_count - D.a_count; double i_per = ((double)i_count / (double)D.tot_count)*100; 
                if (i_per >= a){
                    D.P.root = D.delete_neighbors(G,D.P.root); 
                    D.Initialize_ik_nodes(G,D.P.root); 
                }
                D.a_count = 0; D.tot_count = 1; 
            }    
            auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;    
            avg_time+=elapsed.count();  
        }
        avg_time/=ITERATIONS; time[count] = avg_time; count++;
        std::cout << "Avg. Time: " << avg_time << std::endl;
        std::cout << std::endl; 
    }

    std::ofstream myfile; myfile.open("time_vs_per.txt");

    count = 0; 
    for (int a = START_PER; a <= MAX_PER; a+=DEL_PER){
        myfile << a << " " << time[count] << std::endl;
        count++; 
    }

    myfile.close(); 
        
    return 0;
}