#ifndef HGBEES_H
#define HGBEES_H

#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>
#include <chrono>

const int DIM = 3; 
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
class Cell; 

class Cell{  // Cell - data connected to a given cell
    public:
        double prob;
        std::array <double, DIM> vp; 
        std::array <double, DIM> up; 
        std::array <double, DIM> wp; 
        std::array <double, DIM> vm; 
        std::array <double, DIM> um; 
        std::array <double, DIM> wm; 
        std::array <double, DIM> f; 
        std::array <int, DIM> state; 
        std::array<int, DIM> i_keys; 
        std::array<int, DIM> k_keys; 
        int active; 
        double K; 
};

class Grid{ // Properties of the grid
  public:            
    double thresh;        
    std::array<double,DIM> start;
    std::array<double,DIM> std;
    double dt;
    std::array<double, DIM> del;
    std::array<double, DIM> xh;
    int d;
};

class Lorenz3D{       // Lorenz3D Class
  public:            
    int sigma;        
    int b;
    int r;
    int L;
};

class HGBEES{       // HGBEES Class
  public:
    std::unordered_map<int, Cell> P;
    Cell dead; 

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    void RHS_P(Grid G,Lorenz3D Lor);
    void Record_Data(int step, Grid G);
};

/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
int CantorPair(int state[], int n){
    int key; 
    if(n>2){
        int last = state[n-1]; int new_state[n-1];
        for (int i = 0; i < n - 1; i++){
            new_state[i] = state[i];
        }
        int x = CantorPair(new_state, n-1); int y = last;
        key = (0.5)*(x+y)*(x+y+1)+y;
    }else{
        int x = state[0]; int y = state[1];
        key = (0.5)*(x+y)*(x+y+1)+y;
    }
    return key; 
};

int state_conversion(std::array<int,DIM> state){
    int shift_state[DIM];
    for (int i = 0; i < DIM; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    if(DIM == 1){
        int key = shift_state[0];
        return key;
    }else{
        int key = CantorPair(shift_state, DIM);
        return key;
    }
};

double MC(double th){
    double phi;
    phi = std::max(0.0,std::min({(1+th)/2,2.0,2*th}));
    return phi;
};
/*==============================================================================
Member Function DEFINITIONS
==============================================================================*/
void HGBEES::Initialize_D(Grid G,Lorenz3D Lor){
    Cell blank_c = {.prob = 0, .vp = {0}, .up = {0}, .wp = {0}, .vm = {0}, .um = {0}, .wm = {0}, .f = {0}, .active = -1, .K = 0};
    P[-1] = blank_c; dead = blank_c; 

    std::array<int,DIM> current_state; int current_key; 
    for (int i = round((G.start[0]-2*G.std[0])/G.del[0]); i <= round((G.start[0]+2*G.std[0])/G.del[0]); i++){
        for (int j = round((G.start[1]-2*G.std[1])/G.del[1]); j <= round((G.start[1]+2*G.std[1])/G.del[1]); j++){
            for (int k = round((G.start[2]-2*G.std[2])/G.del[2]); k <= round((G.start[2]+2*G.std[2])/G.del[2]); k++){
                current_state = {i,j,k}; current_key = state_conversion(current_state);

                double x = 0; 
                for(int q = 0; q < DIM; q++){
                    x += pow((current_state[q]*G.del[q]-G.start[q]),2)/pow(G.std[q],2); 
                }
                
                Cell c = {.prob = exp(-4*x/2), .vp = {0}, .up = {0}, .wp = {0}, .vm = {0}, .um = {0}, .wm = {0}, .f = {0}, .state = current_state, .active = 0, .K = 0};
                P[current_key] = c; 
            }
        }
    }
    Initialize_vuw(G,Lor);
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    int current_key; double x[G.d] = {0}; 
    
    for (auto it = P.begin(); it != P.end(); it++){
        current_key = it->first; Cell l_cell = P[current_key]; 
        if(l_cell.active==0){
            std::array<double,DIM> x; 
            for(int i = 0; i < DIM; i++){
                x[i] = G.del[i]*l_cell.state[i];
            }
            
            double v1p = Lor.sigma*(x[1]-(x[0]+G.xh[0]));
            double v2p = -(x[1]+G.xh[1])-x[0]*x[2];
            double v3p = -Lor.b*(x[2]+G.xh[2])+x[0]*x[1]-Lor.b*Lor.r; 
            l_cell.vp = {v1p,v2p,v3p};
            l_cell.up = {std::min(v1p,0.0),std::min(v2p,0.0),std::min(v3p,0.0)};
            l_cell.wp = {std::max(v1p,0.0),std::max(v2p,0.0),std::max(v3p,0.0)}; 
            double v1m = Lor.sigma*(x[1]-(x[0]-G.xh[0]));
            double v2m = -(x[1]-G.xh[1])-x[0]*x[2];
            double v3m = -Lor.b*(x[2]-G.xh[2])+x[0]*x[1]-Lor.b*Lor.r; 
            l_cell.vm = {v1p,v2p,v3p};
            l_cell.um = {std::min(v1m,0.0),std::min(v2m,0.0),std::min(v3m,0.0)};
            l_cell.wm = {std::max(v1m,0.0),std::max(v2m,0.0),std::max(v3m,0.0)}; 
            l_cell.active = 1;
            P[current_key] = l_cell; 
        }
    }
};

void HGBEES::Modify_pointset(Grid G,Lorenz3D Lor){
    double prob_sum = 0; int current_key; double current_prob;  
    int n = P.size(); int keys[n]; int count = 0; 
    for (auto it = P.begin(); it != P.end(); it++) {
        keys[count] = it->first; 
        count++;
    }

    for(int i = 0; i < n; i++){   // Check/Create Neighbors of Big Cells
        current_key = keys[i]; current_prob = P[current_key].prob;
        if(current_prob >= G.thresh){
            prob_sum += P[current_key].prob; 
            std::array<double,DIM> current_v = P[current_key].vp; std::array<int,DIM> current_state = P[current_key].state;    
            
            std::array<int,DIM> num_n; std::array<std::array<int,DIM>,DIM> which_n;  

            for(int c = 0; c < DIM; c++){
                if(current_v[c] > 0){
                    num_n[c] = 2; 
                    which_n[c][0] = current_state[c]; which_n[c][1] = current_state[c]+1;
                }else if(current_v[c]==0){
                    num_n[c] = 3; 
                    which_n[c][0] = current_state[c]-1; which_n[c][1] = current_state[c]; which_n[c][2] = current_state[c]+1;
                }else{
                    num_n[c] = 2;
                    which_n[c][0] = current_state[c]-1; which_n[c][1] = current_state[c];
                }
            }

            for (int i = 0; i < num_n[0]; i++){
                for (int j = 0; j < num_n[1]; j++){
                    for (int k = 0; k < num_n[2]; k++){
                        std::array<int,DIM> new_state = {which_n[0][i], which_n[1][j], which_n[2][k]}; int new_key = state_conversion(new_state);
                        if(!P.count(new_key)){
                            Cell c = {.prob = 0, .vp = {0}, .up = {0}, .wp = {0}, .vm = {0}, .um = {0}, .wm = {0}, .f = {0}, .state = new_state, .active = 0, .K = 0};
                            P[new_key] = c; 
                        }                       
                    }
                }
            }       
        }
    }
    for (int l = 0; l < n; l++){     // Remove small cells that do not neighbor big cells
        int current_key = keys[l];
        if ((P[current_key].prob < G.thresh)&&(P[current_key].active == 1)){
            bool neighbors = true; std::array<double,DIM> current_v = P[current_key].vm; std::array<int,DIM> current_state = P[current_key].state;    

            std::array<int,DIM> num_n; std::array<std::array<int,DIM>,DIM> which_n;  

            for(int c = 0; c < DIM; c++){
                if(current_v[c] < 0){
                    num_n[c] = 2; 
                    which_n[c][0] = current_state[c]; which_n[c][1] = current_state[c]+1;
                }else if(current_v[c]==0){
                    num_n[c] = 3; 
                    which_n[c][0] = current_state[c]-1; which_n[c][1] = current_state[c]; which_n[c][2] = current_state[c]+1;
                }else{
                    num_n[c] = 2;
                    which_n[c][0] = current_state[c]-1; which_n[c][1] = current_state[c];
                }
            }

            for (int i = 0; i < num_n[0]; i++){
                for (int j = 0; j < num_n[1]; j++){
                    for (int k = 0; k < num_n[2]; k++){
                        std::array<int,DIM> new_state = {which_n[0][i], which_n[1][j], which_n[2][k]}; int new_key = state_conversion(new_state);
                        if(P.count(new_key) == 1){
                            if(P[new_key].prob >= G.thresh){
                                neighbors = false;
                            }
                        }                  
                    }
                }
            }   
            if(neighbors){
                P.erase(current_key);
            }
        }
    }

    Initialize_vuw(G,Lor);

    for (auto it = P.begin(); it != P.end(); it++){
        current_key = it->first; 
        P[current_key].prob /= prob_sum;
    } 
};

void HGBEES::RHS_P(Grid G,Lorenz3D Lor){
    for(auto it = P.begin(); it != P.end(); it++){ //Calculating Initial f
        int l_key = it->first; 
        if(l_key != -1){
            Cell l_cell = P[l_key]; std::array<int, DIM> l_state = l_cell.state; l_cell.f = {0.0}; l_cell.K = 0.0;
            for(int q = 0; q < DIM; q++){
                //Initializing i, k nodes
                std::array<int,DIM> i_state = l_state; i_state[q] = i_state[q]-1; int i_key = state_conversion(i_state);
                if(!P.count(i_key)) i_key = -1; l_cell.i_keys[q] = i_key; 
                         
                std::array<int,DIM> k_state = l_state; k_state[q] = k_state[q]+1; int k_key = state_conversion(k_state);
                if(!P.count(k_key)) k_key = -1; Cell k_cell = P[k_key]; l_cell.k_keys[q] = k_key; 

                l_cell.f[q] = l_cell.wp[q] * l_cell.prob + l_cell.up[q] * k_cell.prob;
            }
            P[l_key] = l_cell; 
        }
    }

    for(auto it = P.begin(); it != P.end(); it++){ //Calculating Initial f
        int l_key = it->first;
        if(l_key != -1){
            Cell l_cell = P[l_key]; 
            for(int q = 0; q < DIM; q++){
                int i_key = l_cell.i_keys[q]; Cell i_cell = P[i_key]; 
                if ((l_cell.prob >= G.thresh)||(i_cell.prob >= G.thresh)){ 
                    double F = G.dt*(l_cell.prob-i_cell.prob)/(2*G.del[q]); 
                    for(int e = 0; e < DIM; e++){
                        int j_key; int p_key; Cell j_cell; Cell p_cell; 
                        if (e!=q){
                            j_key = l_cell.i_keys[e]; j_cell = P[j_key];   
                            if (i_key != -1){
                                p_key = i_cell.i_keys[e]; 
                            }else{
                                p_key = -1; 
                            }
                            p_cell = P[p_key]; 
                            
                            l_cell.f[e] -= l_cell.wp[e] * i_cell.wp[q] * F;
                            i_cell.f[e] -= i_cell.wp[e] * i_cell.up[q] * F;
                            j_cell.f[e] -= j_cell.up[e] * i_cell.wp[q] * F;
                            p_cell.f[e] -= p_cell.up[e] * i_cell.up[q] * F; 
                        } 
                        P[j_key] = j_cell; P[p_key] = p_cell;                       
                    }
                    
                    double th,t; 
                    if (i_cell.vp[q]>0){
                        int i_i_key = i_cell.i_keys[q]; Cell i_i_cell = P[i_i_key]; 
                        th = (i_cell.prob-i_i_cell.prob)/(l_cell.prob-i_cell.prob);
                        t = i_cell.vp[q];
                    }else{
                        int k_key = l_cell.k_keys[q]; Cell k_cell = P[k_key]; 
                        th = (k_cell.prob-l_cell.prob)/(l_cell.prob-i_cell.prob);
                        t = -i_cell.vp[q];
                    }
                    
                    i_cell.f[q] += t*(G.del[q]/G.dt - t)*F*MC(th); 
                    P[i_key] = i_cell; 
                }
            }
            P[l_key] = l_cell; 
        }
    }

    for(auto it = P.begin(); it != P.end(); it++){ //Calculating Probability Update
        int l_key = it->first;
        if(l_key != -1){
            Cell l_cell = P[l_key]; std::array<int,DIM> l_state = l_cell.state; 
            for(int q = 0; q < DIM; q++){
                int i_key = l_cell.i_keys[q]; Cell i_cell = P[i_key]; 
                l_cell.K -= (l_cell.f[q]-i_cell.f[q])/G.del[q];  
            }
            l_cell.prob += G.dt*l_cell.K;
            P[l_key] = l_cell;  
        }
    }
};

void HGBEES::Record_Data(int step, Grid G){
    std::string step_str = std::to_string(step); 
    std::string file_name = "pdf_" + step_str + ".txt";
    std::ofstream myfile; myfile.open(file_name);
    for(auto it = P.begin(); it != P.end(); it++){
        int key = it->first; Cell l_cell = P[key]; 
        if (l_cell.prob >= G.thresh){
            myfile << l_cell.prob << " " << l_cell.state[0] << " " << l_cell.state[1] << " " << l_cell.state[2] << std::endl;
        }
    }  
    myfile.close();
};

#endif // HGBEES_H