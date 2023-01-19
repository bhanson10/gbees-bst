#ifndef HGBEES_H
#define HGBEES_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
struct Cell{
    double prob;
    vector<double> v; 
    vector<double> u; 
    vector<double> w; 
    vector<double> f; 
    vector<int> state; 
    int active; 
};

struct Neighbors{
    vector<int> keys; 
    vector<vector<int>> states; 
};

class Grid{       // Grid Class
  public:            
    double thresh;        
    vector<double> start;
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
    vector<int> keys;
    unordered_map<int, Neighbors> Upwind_Neighbors; 
    unordered_map<int, Neighbors> All_Neighbors; 
    int n;        

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    vector<double> RHS_P(Grid G,Lorenz3D Lor);
    void get_keys();
    bool no_neighbors(Grid G, Lorenz3D Lor, int l);
    void Record_Data(string file_name, Grid G);
};

/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
vector<int> ShiftState(vector<int> state, int d){
    vector<int> shift_state(d,0);
    for (int i = 0; i < d; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    return shift_state;
};

int CantorPair(vector<int> state){
    int key; 
    if(state.size()>2){
        int last = state.back(); state.pop_back();
        int x = CantorPair(state); int y = last;
        key = (0.5)*(x+y)*(x+y+1)+y;
    }else{
        int x = state[0]; int y = state[1];
        key = (0.5)*(x+y)*(x+y+1)+y;
    }
    return key; 
};

int state_conversion(vector<int> state, Grid G){
    vector<int> shift_state = ShiftState(state, G.d);
    int key = CantorPair(shift_state);
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
    vector<double> zeros(G.d,0);
    P[-1].prob = 0; P[-1].f = zeros; P[-1].v = zeros; P[-1].u = zeros; P[-1].w = zeros; P[-1].active = 0;
    
    vector<int> current_state;

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                current_state = {i,j,k}; int key = state_conversion(current_state,G);
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[key].prob = exp(-4*x/2); P[key].state = current_state;
            }
        }
    }
    Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    int current_key; vector<int> current_state; vector<double> x(G.d,0); get_keys();
    
    for (int l = 0; l < n; l++){
        current_key = keys[l]; 
        if((!P[current_key].active)&&(current_key!=-1)){
            current_state = P[current_key].state;
            for(int i = 0; i < G.d; i++){
                x[i] = G.del*current_state[i];
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
    double prob_sum = 0; vector<int> dead_keys; 

    for(int l = 0; l < n; l++){   // Check/Create Neighbors of Big Cells
        if (P[keys[l]].prob >= G.thresh){
            int current_key = keys[l]; prob_sum += P[current_key].prob;

            if(Upwind_Neighbors.count(current_key)==0){

                vector<int> current_state = P[current_key].state; vector<vector<int>> n_states; vector<int> n_keys; 

                vector<double> current_v = P[current_key].v;

                vector<int> x_neighbors, y_neighbors, z_neighbors;

                if(current_v[0] > 0) x_neighbors = {current_state[0], current_state[0]+1}; else if (current_v[0] == 0) x_neighbors = {current_state[0]-1,current_state[0],current_state[0]+1}; else x_neighbors = {current_state[0]-1, current_state[0]};
                if(current_v[1] > 0) y_neighbors = {current_state[1], current_state[1]+1}; else if (current_v[1] == 0) y_neighbors = {current_state[1]-1,current_state[1],current_state[1]+1}; else y_neighbors = {current_state[1]-1, current_state[1]};
                if(current_v[2] > 0) z_neighbors = {current_state[2], current_state[2]+1}; else if (current_v[2] == 0) z_neighbors = {current_state[2]-1,current_state[2],current_state[2]+1}; else z_neighbors = {current_state[2]-1, current_state[2]};

                for (int i = 0; i < y_neighbors.size(); i++){
                    for (int j = 0; j < x_neighbors.size(); j++){
                        for (int k = 0; k < z_neighbors.size(); k++){
                            vector<int> new_state = {x_neighbors[j], y_neighbors[i], z_neighbors[k]}; int new_key = state_conversion(new_state,G);
                            n_states.push_back(new_state);
                            n_keys.push_back(new_key);
                            if(!P.count(new_key)){
                                P[new_key].prob = 0; P[new_key].active = 0; P[new_key].state = new_state;
                            }                       
                        }
                    }
                }

                Upwind_Neighbors[current_key].keys = n_keys;  
                Upwind_Neighbors[current_key].states = n_states;  

            }else{
                for (int i = 0; i < Upwind_Neighbors[current_key].keys.size(); i++){
                    int new_key = Upwind_Neighbors[current_key].keys[i]; 

                    if(!P.count(new_key)){
                        P[new_key].prob = 0; P[new_key].active = 0; P[new_key].state = Upwind_Neighbors[current_key].states[i];
                    }
                }
            }
        }else{
            dead_keys.push_back(keys[l]);
        }
    }   

    for (int l = 0; l < dead_keys.size(); l++){     // Remove small cells that do not neighbor big cells
        int current_key = dead_keys[l];
        if(current_key!=-1){
            if (no_neighbors(G,Lor,current_key)){
                P.erase(current_key);
            }
        }
    }
    Initialize_vuw(G,Lor);

    for (int i = 0; i < n; i++){
        P[keys[i]].prob = P[keys[i]].prob/prob_sum;
    }
};

vector<double> HGBEES::RHS_P(Grid G,Lorenz3D Lor){
    vector<double> zeros(G.d,0); vector<double> K(n,0.0); 

    for(int l = 0; l < n; l++){ //Calculating Initial f
        int l_key = keys[l];
        if(l_key != -1){
            vector<int> l_state = P[l_key].state; 
            P[l_key].f = zeros; 
            for(int q = 0; q < G.d; q++){
                vector<int> k_state = l_state; k_state[q] = k_state[q]+1; int k_key = state_conversion(k_state,G);
                if(P.count(k_key)==0) k_key = -1; 
                P[l_key].f[q] = P[l_key].w[q] * P[l_key].prob + P[l_key].u[q] * P[k_key].prob;
            }
        }
    }

    for(int q = 0; q < G.d; q++){
        for(int l = 0; l < n; l++){ // Calculating Total f
            int l_key = keys[l]; 
            if(l_key != -1){
                vector<int> l_state = P[l_key].state; 
                vector<int> i_state = l_state; i_state[q] = i_state[q]-1; 
                int i_key = state_conversion(i_state,G);
                if(P.count(i_key)==0) i_key = -1; 
                if ((P[l_key].prob>G.thresh)||(P[i_key].prob>G.thresh)){
                    double F = G.dt*(P[l_key].prob-P[i_key].prob)/(2*G.del);
                    for(int e = 0; e < G.d; e++){
                        if (e!=q){
                            vector<int> j_state = l_state; j_state[e]=j_state[e]-1; 
                            int j_key = state_conversion(j_state,G); if(P.count(j_key)==0) j_key = -1; 
                            vector<int> p_state = i_state; p_state[e] = p_state[e]-1; 
                            int p_key = state_conversion(p_state,G); if(P.count(p_key)==0) p_key = -1;
                            
                            P[l_key].f[e] -= P[l_key].w[e] * P[i_key].w[q] * F;
                            P[j_key].f[e] -= P[j_key].u[e] * P[i_key].w[q] * F;
                            P[i_key].f[e] -= P[i_key].w[e] * P[i_key].u[q] * F;
                            P[p_key].f[e] -= P[p_key].u[e] * P[i_key].u[q] * F; 
                        }                       
                    }

                    vector<int> i_i_state = i_state; i_i_state[q] = i_i_state[q]-1;
                    int i_i_key = state_conversion(i_i_state,G); if(P.count(i_i_key)==0) i_i_key = -1; 
                    vector<int> k_state = l_state; k_state[q] = k_state[q] + 1; 
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
    for(int l = 0; l < n; l++){
        int l_key = keys[l]; 
        if(l_key != -1){
            vector<int> l_state = P[l_key].state;

            for(int q = 0; q < G.d; q++){
                vector<int> i_state = l_state; i_state[q] = i_state[q] - 1;
                int i_key = state_conversion(i_state,G); if(P.count(i_key)==0) i_key = -1;  
                K[l] -= (P[l_key].f[q]-P[i_key].f[q])/G.del;  
            }
        }
    }
    
    return K;
};

void HGBEES::get_keys(){
    n = P.size(); vector<int> all_keys(n,0); keys = all_keys;

    int count = 0; 
    for (auto it = P.begin(); it != P.end(); it++) {
        keys[count] = it->first;
        count++;
    }
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, int current_key){
    bool neighbors = true; vector<int> current_state = P[current_key].state;
    
    if(All_Neighbors.count(current_key) == 0){
        vector<vector<int>> n_states; vector<int> n_keys; 

        vector<int> x_neighbors = {current_state[0]-1, current_state[0], current_state[0]+1};
        vector<int> y_neighbors = {current_state[1]-1, current_state[1], current_state[1]+1};
        vector<int> z_neighbors = {current_state[2]-1, current_state[2], current_state[2]+1};
        
        for (int i = 0; i < x_neighbors.size(); i++){
            for (int j = 0; j < y_neighbors.size(); j++){
                for (int k = 0; k < z_neighbors.size(); k++){
                    vector<int> new_state = {x_neighbors[i], y_neighbors[j], z_neighbors[k]}; int new_key = state_conversion(new_state,G);
                    n_states.push_back(new_state); n_keys.push_back(new_key); 

                    if(P.count(new_key) == 1){
                        if(P[new_key].prob >= G.thresh){
                            neighbors = false;
                        }
                    }      
                }
            }
        }

        All_Neighbors[current_key].keys = n_keys; 
        All_Neighbors[current_key].states = n_states; 

        return neighbors;

    }else{
        for (int i = 0; i < All_Neighbors[current_key].keys.size(); i++){
            int new_key = All_Neighbors[current_key].keys[i];

            if(P.count(new_key) == 1){
                if(P[new_key].prob >= G.thresh){
                    return neighbors = false;
                }
            }  
        }
    }

    return neighbors; 
};

void HGBEES::Record_Data(string file_name, Grid G){
	ofstream myfile; myfile.open(file_name);
    for(int i = 0; i < n; i++){
        int key = keys[i]; 
        if (P[key].prob >= G.thresh){
            vector<int> current_state = P[key].state;
            myfile << P[key].prob << " " << current_state[0] << " " << current_state[1] << " " << current_state[2] << endl;
        }
    }  
    myfile.close();
};


#endif // HGBEES_H