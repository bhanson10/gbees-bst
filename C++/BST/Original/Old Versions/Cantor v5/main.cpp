/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "GBEES.h"

int main() {

    std::cout << "Initializing Grid and Trajectory..." << std::endl; 
    //=========================== Begin User Input =============================
    int T = 1; Grid G; G.thresh = 0.00002; G.dt = 0.0005; G.del = 0.4;
    double start[] = {-11.5, -10, 9.5}; G.start = start; G.xh = G.del/2; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; Lor.L = 30;
    //============================ End User Input ==============================
    GBEES D;

    std::cout << "Initializing Distribution..." << std::endl; 
    
    D.Initialize_D(G,Lor); D.Modify_pointset(G, Lor); D.Record_Data("pdf_0.txt", G, D.P.root);
    
    std::cout << "Entering Time-marching scheme..." << std::endl; 
    
    for (int i = 1; i <= 5; i++){
        std::cout << D.P.size(D.P.root) << std::endl; 
        D.Modify_pointset(G,Lor); D.RHS_P(G, Lor); D.update_prob(G,D.P.root);
        
        if (i % 400 == 0){ 
            std::cout << "Timestep: " << i << std::endl; 
            std::string filename = "pdf_" + std::to_string(i) + ".txt"; D.Record_Data(filename, G, D.P.root);
        }
    }    

    return 0;
}