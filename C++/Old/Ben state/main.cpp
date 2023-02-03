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
    
    double t = 0; vector<double> K; D.Modify_pointset(G,Lor); D.Record_Data("pdf_0.txt", G);
    
    cout << "Entering Time-marching scheme..." << endl; 

    for (int i = 1; i <= 400; i++){
        t += G.dt; D.Modify_pointset(G,Lor); K = D.RHS_P(G, Lor); 

        for (int j = 0; j < D.n; j++){
            D.P[D.keys[j]].prob = D.P[D.keys[j]].prob + G.dt*K[j];
        }  

        if (i % 400 == 0){ 
            cout << "Timestep: " << i << endl; 
            string filename = "pdf_" + to_string(i) + ".txt"; D.Record_Data(filename, G);
        }
    }    

    return 0;
}