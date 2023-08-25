#ifndef HGBEES_H
#define HGBEES_H

#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

const int DIM = 3; 

/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
struct Cell{
    double prob;
    double* v = new double[DIM]; 
    double* u = new double[DIM];
    double* w = new double[DIM];
    double* f = new double[DIM];
    int* state = new int[DIM];
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
    unordered_map<int, Cell> P;     
    int* keys; 
    int size; 

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    void RHS_P(Grid G,Lorenz3D Lor);
    void get_keys(); 
    bool no_neighbors(Grid G, Lorenz3D Lor, int l);
    void Record_Data(string file_name, Grid G);
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

int state_conversion(int state[], Grid G){
    int shift_state[G.d];
    for (int i = 0; i < G.d; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    int key = CantorPair(shift_state, G.d);
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
    double zeros[G.d] = {0.0}; 
    P[-1].prob = 0; P[-1].f = zeros; P[-1].v = zeros; P[-1].u = zeros; P[-1].w = zeros; P[-1].active = 0;

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                int current_state[] = {i,j,k}; int key = state_conversion(current_state,G); 
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[key].prob = exp(-4*x/2); P[key].active = 0; 
                P[key].state[0] = i; P[key].state[1] = j; P[key].state[2] = k;  
            }
        }
    }
    Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    int current_key; int x[G.d] = {0}; get_keys();

    for (int i = 0; i < size; i++){
        current_key = keys[i];
        if((P[current_key].active==0)&&(current_key!=-1)){
            int* current_state = P[current_key].state;
            for(int i = 0; i < G.d; i++){
                x[i] = G.del*current_state[i];
            }
            double v1 = Lor.sigma*(x[1]-(x[0]+G.xh));
            double v2 = -(x[1]+G.xh)-x[0]*x[2];
            double v3 = -Lor.b*(x[2]+G.xh)+x[0]*x[1]-Lor.b*Lor.r; 
        
            P[current_key].v[0] = v1; P[current_key].v[1] = v2; P[current_key].v[2] = v3; 
            P[current_key].u[0] = min(v1,0.0); P[current_key].u[1] = min(v2,0.0); P[current_key].u[2] = min(v3,0.0); 
            P[current_key].w[0] = max(v1,0.0); P[current_key].w[1] = max(v2,0.0); P[current_key].w[2] = max(v3,0.0); 
            P[current_key].active = 1;
        }
    }
};

void HGBEES::Modify_pointset(Grid G,Lorenz3D Lor){
    double prob_sum = 0; int current_key; 
    
    for(int i = 0; i < size; i++){
        current_key = keys[i];  
        if (P[current_key].prob >= G.thresh){
            prob_sum += P[current_key].prob; 

            int* current_state = P[current_key].state; double* current_v = P[current_key].v;
            
            int x_count, y_count, z_count = 0; int* x_n; int* y_n; int* z_n; 
            if(current_v[0] > 0){
                int x_neighbors[] = {current_state[0], current_state[0]+1};
                x_n = x_neighbors; x_count = 2; 
            }else if (current_v[0] < 0){
                int x_neighbors[] = {current_state[0]-1, current_state[0]};
                x_n = x_neighbors; x_count = 2;
            }else{
                int x_neighbors[] = {current_state[0]-1, current_state[0], current_state[0]+1};
                x_n = x_neighbors; x_count = 3;
            }
            if(current_v[1] > 0){
                int y_neighbors[] = {current_state[1], current_state[1]+1};
                y_n = y_neighbors; y_count = 2;
            }else if (current_v[0] < 0){
                int y_neighbors[] = {current_state[1]-1, current_state[1]};
                y_n = y_neighbors; y_count = 2;
            }else{
                int y_neighbors[] = {current_state[1]-1, current_state[1], current_state[1]+1};
                y_n = y_neighbors; y_count = 3;
            }
            if(current_v[2] > 0){
                int z_neighbors[] = {current_state[2], current_state[2]+1};
                z_n = z_neighbors; z_count = 2;
            }else if (current_v[0] < 0){
                int z_neighbors[] = {current_state[2]-1, current_state[2]};
                z_n = z_neighbors; z_count = 2;
            }else{
                int z_neighbors[] = {current_state[2]-1, current_state[2], current_state[2]+1};
                z_n = z_neighbors; z_count = 3;
            }


            for (int i=0; i < x_count; i++){
                for (int j=0; j < y_count; j++){
                    for (int k=0; k < z_count; k++){
                        int new_state[] = {x_n[i], y_n[j], z_n[k]}; int new_key = state_conversion(new_state,G);
                        if(!P.count(new_key)){
                            P[new_key].prob = 0; P[new_key].active = 0; 
                            P[new_key].state[0] = x_n[i]; P[new_key].state[1] = y_n[j]; P[new_key].state[2] = z_n[k]; 
                        }                       
                    }
                }
            }
        }
    }   
    
    for (int i = 0; i < size; i++){ // Remove small cells that do not neighbor big cells
        current_key = keys[i]; 
        if(current_key!=-1){
            if((P[current_key].prob < G.thresh)&&(no_neighbors(G,Lor,current_key))){
                P.erase(current_key); 
            }
        }
    }
    
    Initialize_vuw(G,Lor);

    for (int i = 0; i < size; i++){
        current_key = keys[i];
        P[current_key].prob = P[current_key].prob/prob_sum;
    }
};

void HGBEES::RHS_P(Grid G,Lorenz3D Lor){
    double zeros[G.d] = {0.0}; 

    for (int l = 0; l < size; l++){
        int l_key = keys[l]; P[l_key].K = 0.0; // Reseting the K value
        if(l_key != -1){
            int* l_state = P[l_key].state; 
            P[l_key].f = zeros; 
            for(int q = 0; q < G.d; q++){
                int* k_state = l_state; k_state[q] = k_state[q]+1; int k_key = state_conversion(k_state,G);
                if(P.count(k_key)==0) k_key = -1; 
                P[l_key].f[q] = P[l_key].w[q] * P[l_key].prob + P[l_key].u[q] * P[k_key].prob;
            }
        }
    }

    for(int q = 0; q < G.d; q++){
        for (int l = 0; l < size; l++){ // Calculating Total f
            int l_key = keys[l]; 
            if(l_key != -1){
                int* l_state = P[l_key].state; 
                int* i_state = l_state; i_state[q] = i_state[q]-1; 
                int i_key = state_conversion(i_state,G);
                if(P.count(i_key)==0) i_key = -1; 
                if ((P[l_key].prob>G.thresh)||(P[i_key].prob>G.thresh)){
                    double F = G.dt*(P[l_key].prob-P[i_key].prob)/(2*G.del);
                    for(int e = 0; e < G.d; e++){
                        if (e!=q){
                            int* j_state = l_state; j_state[e]=j_state[e]-1; 
                            int j_key = state_conversion(j_state,G); if(P.count(j_key)==0) j_key = -1; 
                            int* p_state = i_state; p_state[e] = p_state[e]-1; 
                            int p_key = state_conversion(p_state,G); if(P.count(p_key)==0) p_key = -1;
                            
                            P[l_key].f[e] -= P[l_key].w[e] * P[i_key].w[q] * F;
                            P[j_key].f[e] -= P[j_key].u[e] * P[i_key].w[q] * F;
                            P[i_key].f[e] -= P[i_key].w[e] * P[i_key].u[q] * F;
                            P[p_key].f[e] -= P[p_key].u[e] * P[i_key].u[q] * F; 
                        }                       
                    }

                    int* i_i_state = i_state; i_i_state[q] = i_i_state[q]-1;
                    int i_i_key = state_conversion(i_i_state,G); if(P.count(i_i_key)==0) i_i_key = -1; 
                    int* k_state = l_state; k_state[q] = k_state[q] + 1; 
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
    
    for (int l = 0; l < size; l++){ 
        int l_key = keys[l]; 
        if(l_key != -1){
            int* l_state = P[l_key].state;

            for(int q = 0; q < G.d; q++){
                int* i_state = l_state; i_state[q] = i_state[q] - 1;
                int i_key = state_conversion(i_state,G); if(P.count(i_key)==0) i_key = -1;  
                P[l_key].K -= (P[l_key].f[q]-P[i_key].f[q])/G.del;  
            }
        }
    }
};

void HGBEES::get_keys(){
    size = P.size();  int all_keys[size] = {0}; int count = 0; 
    
    for (auto& it: P){
        all_keys[count] = it.first;
        count++;
    }
    keys = all_keys;
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, int current_key){
    bool neighbors = true; int* current_state = P[current_key].state;

    int x_neighbors[] = {current_state[0]-1, current_state[0], current_state[0]+1};
    int y_neighbors[] = {current_state[1]-1, current_state[1], current_state[1]+1};
    int z_neighbors[] = {current_state[2]-1, current_state[2], current_state[2]+1};
    
    for (int i = 0; i < G.d; i++){
        for (int j = 0; j < G.d; j++){
            for (int k = 0; k < G.d; k++){
                int new_state[] = {x_neighbors[i], y_neighbors[j], z_neighbors[k]}; int new_key = state_conversion(new_state,G);

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
    for (auto& it: P){ 
        int key = it.first;
        if (P[key].prob >= G.thresh){
            int* current_state = P[key].state;
            myfile << P[key].prob << " " << current_state[0] << " " << current_state[1] << " " << current_state[2] << endl;
        }
    }  
    myfile.close();
};


#endif // HGBEES_H