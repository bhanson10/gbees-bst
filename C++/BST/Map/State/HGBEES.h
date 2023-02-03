#ifndef HGBEES_H
#define HGBEES_H

#include <iostream>
#include <map>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
using namespace std;

const int DIM = 3; // Dimension
const int MAX = 50000; // Max dictionary size 
array<int,DIM> DEAD = {-1,-1,-1}; //Global Dead State
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
struct Cell{  // Cell Class
    double prob;
    array<double, DIM> v;
    array<double, DIM> u;
    array<double, DIM> w;
    array<double, DIM> f;
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
    map<array<int,DIM>, Cell> P;

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    void RHS_P(Grid G,Lorenz3D Lor);
    bool no_neighbors(Grid G, Lorenz3D Lor, array<int,DIM> l);
    void Record_Data(string file_name, Grid G);
};

/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
double MC(double th){
    double phi;
    phi = max(0.0,min({(1+th)/2,2.0,2*th}));
    return phi;
};
/*==============================================================================
Member Function DEFINITIONS
==============================================================================*/
void HGBEES::Initialize_D(Grid G,Lorenz3D Lor){
    P[DEAD].active = 0; P[DEAD].v = {0.0}; P[DEAD].u = {0.0}; P[DEAD].w = {0.0}; P[DEAD].f = {0.0};
    
    array<int,DIM> current_state; 

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                current_state = {i,j,k};
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[current_state].prob = exp(-4*x/2); P[current_state].active = 0; 
            }
        }
    }
    Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    array<int,DIM> current_key; double x[G.d] = {0}; 
    
    for (auto it = P.begin(); it != P.end(); it++){
        current_key = it->first; 
        if((P[current_key].active==0)&&(current_key!=DEAD)){
            for(int i = 0; i < G.d; i++){
                x[i] = G.del*current_key[i]; 
            }
            
            double v1 = Lor.sigma*(x[1]-(x[0]+G.xh));
            double v2 = -(x[1]+G.xh)-x[0]*x[2];
            double v3 = -Lor.b*(x[2]+G.xh)+x[0]*x[1]-Lor.b*Lor.r; 
            P[current_key].v = {v1,v2,v3};
            P[current_key].u = {min(v1,0.0),min(v2,0.0),min(v3,0.0)};
            P[current_key].w = {max(v1,0.0),max(v2,0.0),max(v3,0.0)}; 
            P[current_key].active = 1;
        }
    }
};

void HGBEES::Modify_pointset(Grid G,Lorenz3D Lor){
    double prob_sum = 0; array<int,DIM> current_state; double current_prob; array<double,DIM> current_v; array<int,DIM> new_state; 
    int added = 0; 

    int n = 0; array<array<int,DIM>,50000> all_keys; // Getting keys that will be operated on
    for (auto it = P.begin(); it != P.end(); it++) {
        all_keys[n] = it->first; 
        n++;
    }

    for(int l = 0; l < n; l++){   // Check/Create Neighbors of Big Cells
        current_state = all_keys[l]; current_prob = P[current_state].prob;
        if (current_prob >= G.thresh){
            prob_sum += current_prob; 

            current_v = P[current_state].v;

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
                        new_state = {x_n[i], y_n[j], z_n[k]}; 
                        if(!P.count(new_state)){
                            P[new_state].prob = 0; P[new_state].active = 0; 
                            added++; 
                        }                       
                    }
                }
            }               
        }
    }   
    for (int l = 0; l < n; l++){     // Remove small cells that do not neighbor big cells
        current_state = all_keys[l];
        if(current_state!=DEAD){
            if ((P[current_state].prob < G.thresh)&&(no_neighbors(G,Lor,current_state))){
                P.erase(current_state);
            }
        }
    }
    
    Initialize_vuw(G,Lor);

    for (auto it = P.begin(); it != P.end(); it++){
        current_state = it->first; 
        P[current_state].prob /= prob_sum;
    } 
};


void HGBEES::RHS_P(Grid G,Lorenz3D Lor){
    array<int,DIM> l_key; array<int,DIM> k_key; array<int,DIM> p_key; 
    array<int,DIM> j_key; array<int,DIM> i_i_key; array<int,DIM> i_key; 

    for(auto it = P.begin(); it != P.end(); it++){ //Calculating Initial f
        l_key = it->first;
        if(l_key != DEAD){
            P[l_key].K = 0.0; P[l_key].f = {0.0};
            for(int q = 0; q < G.d; q++){
                k_key = l_key; k_key[q] = k_key[q]+1; 
                if(P.count(k_key)==0) k_key = DEAD; 
                P[l_key].f[q] = P[l_key].w[q] * P[l_key].prob + P[l_key].u[q] * P[k_key].prob;
            }
        }
    }

    for(int q = 0; q < G.d; q++){
        for(auto it = P.begin(); it != P.end(); it++){ // Calculating Total f
            l_key = it->first; 
            if(l_key != DEAD){
                i_key = l_key; i_key[q] = i_key[q]-1; if(P.count(i_key)==0) i_key = DEAD; 
                if ((P[l_key].prob>G.thresh)||(P[i_key].prob>G.thresh)){
                    double F = G.dt*(P[l_key].prob-P[i_key].prob)/(2*G.del);
                    for(int e = 0; e < G.d; e++){
                        if (e!=q){
                            j_key = l_key; j_key[e] -= 1; if(P.count(j_key)==0) j_key = DEAD;
                            p_key = i_key; p_key[e] -= 1; if(P.count(p_key)==0) p_key = DEAD;
                            
                            P[l_key].f[e] -= P[l_key].w[e] * P[i_key].w[q] * F;
                            P[j_key].f[e] -= P[j_key].u[e] * P[i_key].w[q] * F;
                            P[i_key].f[e] -= P[i_key].w[e] * P[i_key].u[q] * F;
                            P[p_key].f[e] -= P[p_key].u[e] * P[i_key].u[q] * F; 
                        }                       
                    }
                    i_i_key = i_key; i_i_key[q]-=1; if(P.count(i_i_key)==0) i_i_key = DEAD; 
                    k_key = l_key; k_key[q] += 1; if(P.count(k_key)==0) k_key = DEAD;
                    
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
        l_key = it->first; 
        if(l_key != DEAD){
            for(int q = 0; q < G.d; q++){
                i_key = l_key; i_key[q] -= 1; if(P.count(i_key)==0) i_key = DEAD;  
                P[l_key].K -= (P[l_key].f[q]-P[i_key].f[q])/G.del;  
            }
        }
    }
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, array<int,DIM> current_state){
    bool neighbors = true; array<int,DIM> new_state; 

    for (int i = current_state[0]-1; i <= current_state[0]+1; i++){
        for (int j = current_state[1]-1; j <= current_state[1]+1; j++){
            for (int k = current_state[2]-1; k <= current_state[2]+1; k++){
                new_state = {i, j, k}; 
                if(P.count(new_state) == 1){
                    if(P[new_state].prob >= G.thresh){
                        return neighbors = false;
                    }
                }                    
            }
        }
    }   
    return neighbors;
};

void HGBEES::Record_Data(string file_name, Grid G){
	ofstream myfile; myfile.open(file_name); array<int,DIM> key; 
    for(auto it = P.begin(); it != P.end(); it++){
        key = it->first; 
        if (P[key].prob >= G.thresh){
            myfile << P[key].prob << " " << key[0] << " " << key[1] << " " << key[2] << endl;
        }
    }  
    myfile.close();
};
#endif // HGBEES_H