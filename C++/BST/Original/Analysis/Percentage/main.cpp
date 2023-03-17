/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main() {
 
    //===================================== Begin User Input ======================================
    Grid G; G.thresh = 0.00002; G.dt = 0.0005; G.del = {0.4,0.4,0.4}; G.start = {-11.5, -10, 9.5};
    G.std = {1.0, 1.0, 1.0}; G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2}; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; Lor.L = 30;
    //====================================== End User Input =======================================
    const int MAX_PER = 75; 
    const int START_PER = 5; 
    const int DEL_PER = 5; 
    const int ITERATIONS = 3; 
    const int GBEES_STEPS = 2000; 
    const int LENGTH = ((MAX_PER-START_PER)/DEL_PER)+1;
    int count = 0; std::array<double, LENGTH> time; 

    for (int a = START_PER; a <= MAX_PER; a+=DEL_PER){
        double avg_time = 0; 
        std::cout << "Percentage: " << a << "%" << std::endl; 
        for (int b = 1; b <= ITERATIONS; b++){

            GBEES D;
            
            D.Initialize_D(G,Lor); D.Modify_pointset(G, Lor); 

            auto start = std::chrono::high_resolution_clock::now(); 
            std::chrono::duration<double> elapsed;

            for (int i = 1; i <= GBEES_STEPS; i++){
                D.Modify_pointset(G,Lor); D.RHS_P(G, Lor);

                if ((1-(D.P.active_size(D.P.root,G)/D.P.full_size(D.P.root)))*100 > a){
                    D.P.root = D.delete_neighbors(G,D.P.root); 
                }
            }
            auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;    
            avg_time+=elapsed.count();  
        }
        avg_time/=ITERATIONS; time[count] = avg_time; count++;
        std::cout << "Avg. Time: " << avg_time << std::endl;
        std::cout << std::endl; 
    }

    std::ofstream myfile; myfile.open("time_vs_percent.txt");

    for (int a = 0; a < LENGTH; a++){
        myfile << (a+1)*DEL_PER << " " << time[a] << std::endl;
    }

    myfile.close(); 

    return 0;
}