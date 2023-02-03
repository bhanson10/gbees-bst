/*==============================================================================

MAIN.CPP

==============================================================================*/
#include "HGBEES.h"

int main() {

    cout << "Initializing Grid and Trajectory..." << endl; 
    //=========================== Begin User Input =============================
    int T = 1; Grid G; G.thresh = 0.00002; G.dt = 0.0005; G.d = 3;
    G.start = {-11.5, -10, 9.5}; G.del = 0.4; G.xh = G.del/2; 
    Lorenz3D Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; Lor.L = 30;
    //============================ End User Input ==============================
    HGBEES D;

    cout << "Initializing Distribution..." << endl; 
    
    D.Initialize_D(G,Lor); 
    
    D.Modify_pointset(G,Lor); D.Record_Data("pdf_0.txt", G);

    cout << "Entering Time-marching scheme..." << endl; 

    for (int i = 1; i <= 400; i++){
        D.Modify_pointset(G,Lor); D.RHS_P(G, Lor); 

        int current_key; 
        for (auto& it: D.P){
            current_key = it.first;
            D.P[current_key].prob += G.dt*D.P[current_key].K;
        }  

        if (i % 400 == 0){ 
            cout << "Timestep: " << i << endl; 
            string filename = "pdf_" + to_string(i) + ".txt"; D.Record_Data(filename, G);
        }
    }    

    return 0;

}