/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "HGBEES.h"

int main() {

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    Grid G; G.thresh = 0.00002; G.dt = 0.0005; G.del = {0.4,0.4,0.4}; G.start = {-11.5, -10, 9.5};
    G.std = {1.0, 1.0, 1.0}; G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2}; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; Lor.L = 30;
    //====================================== End User Input =======================================
    HGBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,Lor); D.Modify_pointset(G, Lor); D.Record_Data(0,G); 
    
    std::cout << "Entering Time-marching scheme..." << std::endl; 
    
    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;
    
    for (int i = 1; i <= 5; i++){
        std::cout << D.P.size() << std::endl; 
        D.Modify_pointset(G,Lor); D.RHS_P(G, Lor);
        
        if (i % 400 == 0){ 
            auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
            std::cout << "Timestep: " << i << ", Time: " << elapsed.count() << " s, Total Cells: " << D.P.size() << std::endl;
            D.Record_Data(i, G);
        }
        
    }    

    return 0;
}