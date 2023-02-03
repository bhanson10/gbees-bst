#ifndef HGBEES_H
#define HGBEES_H

#include <iostream>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
using namespace std;

const int DIM = 3; 
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
struct Cell{  // Cell Class
    double prob;
    double v[DIM]; 
    double u[DIM]; 
    double w[DIM]; 
    double f[DIM]; 
    int state[DIM]; 
    double K; 
    int active; 
};

class Grid{       // Grid Class
  public:            
    double thresh;        
    double* start;
    double dt;
    double del;
    double xh;
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
    map<int, Cell> P;

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    void RHS_P(Grid G,Lorenz3D Lor);
    bool no_neighbors(Grid G, Lorenz3D Lor, int l);
    void Record_Data(string file_name, Grid G);
};

/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
int RosenbergPair(int state[], int d, int m){
    int key; 
    if(d == 1){
        return key = state[0];
    }else{
        int new_state[d-1];
        for (int i = 0; i < d-1; i++){
            new_state[i] = state[i]; 
        }
        int new_m = *max_element(new_state, new_state + d-1);
        return key = RosenbergPair(new_state, d-1, new_m) + pow(m,d) + (m - state[d-1])*(pow(m+1, d-1) - pow(m,d-1)); 
    }
};

int state_conversion(int state[], Grid G){
    int shift_state[G.d];
    for (int i = 0; i < G.d; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    int m = *max_element(shift_state, shift_state + G.d);
    int key = RosenbergPair(shift_state, G.d, m);
    return key;
};

double MC(double th){
    double phi;
    phi = max(0.0,min({(1+th)/2,2.0,2*th}));
    return phi;
};
/*==============================================================================
Member Function DEFINITIONS
==============================================================================*/
void HGBEES::Initialize_D(Grid G,Lorenz3D Lor){
    double zeros[G.d] = {0}; P[-1].active = 0;
    memcpy(P[-1].f, zeros, sizeof(zeros)); memcpy(P[-1].v, zeros, sizeof(zeros));
    memcpy(P[-1].u, zeros, sizeof(zeros)); memcpy(P[-1].w, zeros, sizeof(zeros));

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                int current_state[] = {i,j,k}; int key = state_conversion(current_state,G);
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[key].prob = exp(-4*x/2); P[key].active = 0; 
                memcpy(P[key].state, current_state, sizeof(current_state));  
            }
        }
    }
    Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    int current_key; double x[G.d] = {0}; 
    
    for (auto it = P.begin(); it != P.end(); it++){
        current_key = it->first;
        if((P[current_key].active==0)&&(current_key!=-1)){
            for(int i = 0; i < G.d; i++){
                x[i] = G.del*P[current_key].state[i];
            }
            
            double v1 = Lor.sigma*(x[1]-(x[0]+G.xh));
            double v2 = -(x[1]+G.xh)-x[0]*x[2];
            double v3 = -Lor.b*(x[2]+G.xh)+x[0]*x[1]-Lor.b*Lor.r; 
            double total_v[] = {v1,v2,v3};
            double total_u[] = {min(v1,0.0),min(v2,0.0),min(v3,0.0)};
            double total_w[] = {max(v1,0.0),max(v2,0.0),max(v3,0.0)}; 

            memcpy(P[current_key].v, total_v, sizeof(total_v));
            memcpy(P[current_key].u, total_u, sizeof(total_u));
            memcpy(P[current_key].w, total_w, sizeof(total_w));
            P[current_key].active = 1;
        }
    }
};

void HGBEES::Modify_pointset(Grid G,Lorenz3D Lor){
    double prob_sum = 0; int current_key; double current_prob; double current_v[DIM]; int current_state[DIM];

    int count = 0; int n = P.size(); int all_keys[n]; // Getting keys that will be operated on
    for (auto it = P.begin(); it != P.end(); it++) {
        all_keys[count] = it->first;
        count++;
    }

    for(int l = 0; l < n; l++){   // Check/Create Neighbors of Big Cells
        current_key = all_keys[l]; current_prob = P[current_key].prob;
        if (current_prob >= G.thresh){
            prob_sum += current_prob; 

            memcpy(current_v, P[current_key].v, sizeof(P[current_key].v));
            memcpy(current_state, P[current_key].state, sizeof(P[current_key].state));

            int x_c = 0; int y_c = 0; int z_c = 0; 
            int x_n[DIM]; int y_n[DIM]; int z_n[DIM];

            if(current_v[0] > 0){ //Getting x neighbors
                x_c = 2; x_n[0] = current_state[0]; x_n[1] = current_state[0]+1;
            }else if (current_v[0] == 0){
                x_c = 3; x_n[0] = current_state[0]-1; x_n[1] = current_state[0]; x_n[2] = current_state[0]+1;
            }else{
                x_c = 2; x_n[0] = current_state[0]-1; x_n[1] = current_state[0];
            }

            if(current_v[1] > 0){ //Getting y neighbors
                y_c = 2; y_n[0] = current_state[1]; y_n[1] = current_state[1]+1;
            }else if (current_v[1] == 0){
                y_c = 3; y_n[0] = current_state[1]-1; y_n[1] = current_state[1]; y_n[2] = current_state[1]+1;
            }else{
                y_c = 2; y_n[0] = current_state[1]-1; y_n[1] = current_state[1];
            }

            if(current_v[2] > 0){ //Getting z neighbors
                z_c = 2; z_n[0] = current_state[2]; z_n[1] = current_state[2]+1;
            }else if (current_v[2] == 0){
                z_c = 3; z_n[0] = current_state[2]-1; z_n[1] = current_state[2]; z_n[2] = current_state[2]+1;
            }else{
                z_c = 2; z_n[0] = current_state[2]-1; z_n[1] = current_state[2];
            }

            for (int i = 0; i < x_c; i++){
                for (int j = 0; j < y_c; j++){
                    for (int k = 0; k < z_c; k++){
                        int new_state[] = {x_n[i], y_n[j], z_n[k]}; int new_key = state_conversion(new_state,G);
                        if(!P.count(new_key)){
                            P[new_key].prob = 0; P[new_key].active = 0; 
                            memcpy(P[new_key].state, new_state, sizeof(new_state));
                        }                       
                    }
                }
            }               
        }
    }   

    for (int l = 0; l < n; l++){     // Remove small cells that do not neighbor big cells
        int current_key = all_keys[l];
        if(current_key!=-1){
            if ((P[current_key].prob < G.thresh)&&(no_neighbors(G,Lor,current_key))){
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
    double zeros[G.d] = {0}; int l_state[DIM]; int k_state[DIM]; int p_state[DIM]; int j_state[DIM]; 
    int i_i_state[DIM]; int i_state[DIM]; 

    for(auto it = P.begin(); it != P.end(); it++){ //Calculating Initial f
        int l_key = it->first;
        if(l_key != -1){
            memcpy(l_state, P[l_key].state, sizeof(P[l_key].state)); P[l_key].K = 0.0;
            memcpy(P[l_key].f, zeros, sizeof(zeros)); 
            for(int q = 0; q < G.d; q++){
                memcpy(k_state, l_state, sizeof(l_state)); k_state[q] += 1; int k_key = state_conversion(k_state,G);
                if(P.count(k_key)==0) k_key = -1; 
                P[l_key].f[q] = P[l_key].w[q] * P[l_key].prob + P[l_key].u[q] * P[k_key].prob;
            }
        }
    }
    
    for(int q = 0; q < G.d; q++){
        for(auto it = P.begin(); it != P.end(); it++){ // Calculating Total f
            int l_key = it->first; 
            if(l_key != -1){
                memcpy(l_state, P[l_key].state, sizeof(P[l_key].state)); 
                memcpy(i_state, l_state, sizeof(l_state)); i_state[q] -= 1; 
                int i_key = state_conversion(i_state,G);
                if(P.count(i_key)==0) i_key = -1; 
                if ((P[l_key].prob>G.thresh)||(P[i_key].prob>G.thresh)){
                    double F = G.dt*(P[l_key].prob-P[i_key].prob)/(2*G.del);
                    for(int e = 0; e < G.d; e++){
                        if (e!=q){
                            memcpy(j_state, l_state, sizeof(l_state)); j_state[e]-=1; 
                            int j_key = state_conversion(j_state,G); if(P.count(j_key)==0) j_key = -1; 
                            memcpy(p_state, i_state, sizeof(i_state)); p_state[e]-=1; 
                            int p_key = state_conversion(p_state,G); if(P.count(p_key)==0) p_key = -1;
                            
                            P[l_key].f[e] -= P[l_key].w[e] * P[i_key].w[q] * F;
                            P[j_key].f[e] -= P[j_key].u[e] * P[i_key].w[q] * F;
                            P[i_key].f[e] -= P[i_key].w[e] * P[i_key].u[q] * F;
                            P[p_key].f[e] -= P[p_key].u[e] * P[i_key].u[q] * F; 
                        }                       
                    }

                    memcpy(i_i_state, i_state, sizeof(i_state)); i_i_state[q]-=1;
                    int i_i_key = state_conversion(i_i_state,G); if(P.count(i_i_key)==0) i_i_key = -1; 
                    memcpy(k_state, l_state, sizeof(l_state)); k_state[q]+=1; 
                    int k_key = state_conversion(k_state,G); if(P.count(k_key)==0) k_key = -1; 
                    
                    double th,t; 
                    if (P[i_key].v[q]>0){
                        th = (P[i_key].prob-P[i_i_key].prob)/(P[l_key].prob-P[i_key].prob);
                    }else{
                        th = (P[k_key].prob-P[l_key].prob)/(P[l_key].prob-P[i_key].prob);
                    }

                    if(P[i_key].v[q]>=0){
                        t = P[i_key].v[q];
                    }else{
                        t = -P[i_key].v[q];
                    } 
                    
                    P[i_key].f[q] += t*(G.del/G.dt - t)*F*MC(th); 
                }
            }
        }
    }
    
    for(auto it = P.begin(); it != P.end(); it++){
        int l_key = it->first; 
        if(l_key != -1){
            memcpy(l_state, P[l_key].state, sizeof(P[l_key].state)); 

            for(int q = 0; q < G.d; q++){
                memcpy(i_state, l_state, sizeof(l_state)); i_state[q]-=1;
                int i_key = state_conversion(i_state,G); if(P.count(i_key)==0) i_key = -1;  
                P[l_key].K -= (P[l_key].f[q]-P[i_key].f[q])/G.del;  
            }
        }
    }
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, int current_key){
    bool neighbors = true; int current_state[DIM]; memcpy(current_state, P[current_key].state, sizeof(P[current_key].state));

    for (int i = current_state[0]-1; i <= current_state[0]+1; i++){
        for (int j = current_state[1]-1; j <= current_state[1]+1; j++){
            for (int k = current_state[2]-1; k <= current_state[2]+1; k++){
                int new_state[] = {i, j, k}; int new_key = state_conversion(new_state,G);
                if(P.count(new_key) == 1){
                    if(P[new_key].prob >= G.thresh){
                        return neighbors = false;
                    }
                }                    
            }
        }
    }   
    return neighbors;
};

void HGBEES::Record_Data(string file_name, Grid G){
	ofstream myfile; myfile.open(file_name);
    for(auto it = P.begin(); it != P.end(); it++){
        int key = it->first; 
        if (P[key].prob >= G.thresh){
            myfile << P[key].prob << " " << P[key].state[0] << " " << P[key].state[1] << " " << P[key].state[2] << endl;
        }
    }  
    myfile.close();
};


#endif // HGBEES_H