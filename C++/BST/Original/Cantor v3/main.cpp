/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main() {

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //===================================== Begin User Input ======================================
    Grid G; G.thresh = 0.00002; G.dt = 0.0005; G.start = {-11.5, -10, 9.5}; G.std = {1.0, 1.0, 1.0}; 
    G.del = {G.std[0]/2,G.std[1]/2,G.std[2]/2}; G.xh = {G.del[0]/2, G.del[1]/2, G.del[2]/2}; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48;
    //====================================== End User Input =======================================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,Lor); D.Modify_pointset(G,Lor); D.Record_Data("pdf_0.txt", G, D.P.root); 
    std::cout << "Timestep: 0, Time: 0 s, Active/Total Cells: " << D.P.active_size(D.P.root, G) << "/" << D.P.full_size(D.P.root) << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;

    for (int i = 1; i <= 2000; i++){
        D.Modify_pointset(G,Lor); D.RHS_P(G, Lor);

        if (i % 20 == 0){
            D.P.root = D.delete_neighbors(G,D.P.root); 
        }
        
        if (i % 400 == 0){ 
            auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
            std::cout << "Timestep: " << i << ", Time: " << elapsed.count() << " s, Active/Total Cells: " << D.P.active_size(D.P.root, G) << "/" << D.P.full_size(D.P.root) << std::endl;
            std::string filename = "pdf_" + std::to_string(i) + ".txt"; D.Record_Data(filename, G, D.P.root);
        }       
    }    

    return 0;
}