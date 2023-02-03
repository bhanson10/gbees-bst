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

vector<int> DEAD = {-1,-1,-1}; //Global Dead State

/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
struct Cell{
    double prob;
    vector<double> v; 
    vector<double> u; 
    vector<double> w; 
    vector<double> f; 
    int active; 
};

struct VectorHasher {
    int operator()(const vector<int> &V) const {
        int hash = V.size();
        for(auto &i : V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
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
    unordered_map<vector<int>, Cell, VectorHasher> P;
    vector<vector<int>> keys;
    unordered_map<vector<int>, vector<vector<int>>, VectorHasher> Upwind_Neighbors; 
    unordered_map<vector<int>, vector<vector<int>>, VectorHasher> All_Neighbors; 
    int n;        

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    vector<double> RHS_P(Grid G,Lorenz3D Lor);
    void get_keys();
    bool no_neighbors(Grid G, Lorenz3D Lor, vector<int> key);
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
    vector<double> zeros(G.d,0);
    P[DEAD].prob = 0; P[DEAD].f = zeros; P[DEAD].v = zeros; P[DEAD].u = zeros; P[DEAD].w = zeros; P[DEAD].active = 0;
    
    vector<int> state;

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                state = {i,j,k};
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[state].prob = exp(-4*x/2);
            }
        }
    }
    get_keys(); Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    vector<int> current_key; vector<double> x(G.d,0); n = P.size(); 
    
    for (int l = 0; l < n; l++){
        current_key = keys[l]; 
        if((!P[current_key].active)&&(current_key != DEAD)){

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
    get_keys();  n = P.size(); 

    for(int l = 0; l < n; l++){   // Check/Create Neighbors of Big Cells
        if (P[keys[l]].prob >= G.thresh){
            vector<int> current_key = keys[l];

            if(Upwind_Neighbors.count(current_key)==0){
            
                vector<double> current_v = P[current_key].v; vector<vector<int>> neighbors;

                vector<int> x_neighbors, y_neighbors, z_neighbors;

                if(current_v[0] > 0) x_neighbors = {current_key[0], current_key[0]+1}; else if (current_v[0] == 0) x_neighbors = {current_key[0]-1,current_key[0],current_key[0]+1}; else x_neighbors = {current_key[0]-1, current_key[0]};
                if(current_v[1] > 0) y_neighbors = {current_key[1], current_key[1]+1}; else if (current_v[1] == 0) y_neighbors = {current_key[1]-1,current_key[1],current_key[1]+1}; else y_neighbors = {current_key[1]-1, current_key[1]};
                if(current_v[2] > 0) z_neighbors = {current_key[2], current_key[2]+1}; else if (current_v[2] == 0) z_neighbors = {current_key[2]-1,current_key[2],current_key[2]+1}; else z_neighbors = {current_key[2]-1, current_key[2]};

                for (int i = 0; i < y_neighbors.size(); i++){
                    for (int j = 0; j < x_neighbors.size(); j++){
                        for (int k = 0; k < z_neighbors.size(); k++){
                            vector<int> new_key = {x_neighbors[j], y_neighbors[i], z_neighbors[k]}; 
                            neighbors.push_back(new_key);
                            if(!P.count(new_key)){
                                P[new_key].prob = 0; P[new_key].active = 0; 
                            }                       
                        }
                    }
                } 

                Upwind_Neighbors[current_key] = neighbors;
            }else{
                for (int i = 0; i < Upwind_Neighbors[current_key].size(); i++){
                    vector<int> new_key = Upwind_Neighbors[current_key][i]; 

                    if(!P.count(new_key)){
                        P[new_key].prob = 0; P[new_key].active = 0;
                    }
                }
            }                
        }
    }   
    n = P.size(); get_keys(); 

    for (int l = 0; l < n; l++){     // Remove small cells that do not neighbor big cells in the upwind direction
        vector<int> current_key = keys[l];
        if(current_key != DEAD){
            if ((P[current_key].active==1)&&(P[current_key].prob < G.thresh)&&(no_neighbors(G,Lor,current_key))){
                P.erase(current_key);
            }
        }
    }
    n = P.size(); get_keys(); double prob_sum;
    for (int i = 0; i < n; i++){
        P[keys[i]].prob = max(P[keys[i]].prob,0.); prob_sum += P[keys[i]].prob;
    }
    for (int i = 0; i < n; i++){
        P[keys[i]].prob = P[keys[i]].prob/prob_sum;
    }
    get_keys(); Initialize_vuw(G,Lor);
};

vector<double> HGBEES::RHS_P(Grid G,Lorenz3D Lor){
    n = P.size(); get_keys(); vector<double> zeros(G.d,0); vector<double> K(n,0.0); 

    for(int l = 0; l < n; l++){ //Calculating Initial f
        vector<int> l_key = keys[l];
        if(l_key != DEAD){
            P[l_key].f = zeros; 
            for(int q = 0; q < G.d; q++){
                vector<int> k_key = l_key; k_key[q] = k_key[q]+1;
                if(P.count(k_key)==0) k_key = DEAD; 
                P[l_key].f[q] = P[l_key].w[q] * P[l_key].prob + P[l_key].u[q] * P[k_key].prob;
            }
        }
    }

    for(int q = 0; q < G.d; q++){
        for(int l = 0; l < n; l++){ // Calculating Total f
            vector<int> l_key = keys[l]; 
            if(l_key != DEAD){
                vector<int> i_key = l_key; i_key[q] = i_key[q]-1; if(P.count(i_key)==0) i_key = DEAD; 
                if ((P[l_key].prob>G.thresh)||(P[i_key].prob>G.thresh)){
                    double F = G.dt*(P[l_key].prob-P[i_key].prob)/(2*G.del);
                    for(int e = 0; e < G.d; e++){
                        if (e!=q){
                            vector<int> j_key = l_key; j_key[e] = j_key[e]-1; if(P.count(j_key)==0) j_key = DEAD; 
                            vector<int> p_key = i_key; p_key[e] = p_key[e]-1; if(P.count(p_key)==0) p_key = DEAD;
                            
                            P[l_key].f[e] -= P[l_key].w[e] * P[i_key].w[q] * F;
                            P[j_key].f[e] -= P[j_key].u[e] * P[i_key].w[q] * F;
                            P[i_key].f[e] -= P[i_key].w[e] * P[i_key].u[q] * F;
                            P[p_key].f[e] -= P[p_key].u[e] * P[i_key].u[q] * F; 
                        }                       
                    }

                    vector<int> i_i_key = i_key; i_i_key[q] = i_i_key[q]-1; if(P.count(i_i_key)==0) i_i_key = DEAD; 
                    vector<int> k_key = l_key; k_key[q] = k_key[q] + 1; if(P.count(k_key)==0) k_key = DEAD; 
             
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
        vector<int> l_key = keys[l]; 
        if(l_key != DEAD){
            for(int q = 0; q < G.d; q++){
                vector<int> i_key = l_key; i_key[q] = i_key[q] - 1; if(P.count(i_key)==0) i_key = DEAD;  
                K[l] -= (P[l_key].f[q]-P[i_key].f[q])/G.del;  
            }
        }
    }
    
    return K;
};

void HGBEES::get_keys(){
    keys.clear();

    for (auto it = P.begin(); it != P.end(); it++) {
        keys.push_back(it->first);
    }
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, vector<int> current_key){
    bool no_neighbors = true;

    if(All_Neighbors.count(current_key) == 0){

        vector<vector<int>> n_keys;  

        vector<int> x_neighbors = {current_key[0]-1, current_key[0], current_key[0]+1};
        vector<int> y_neighbors = {current_key[1]-1, current_key[1], current_key[1]+1};
        vector<int> z_neighbors = {current_key[2]-1, current_key[2], current_key[2]+1};
        
        for (int i = 0; i < x_neighbors.size(); i++){
            for (int j = 0; j < y_neighbors.size(); j++){
                for (int k = 0; k < z_neighbors.size(); k++){
                    vector<int> new_key = {x_neighbors[i], y_neighbors[j], z_neighbors[k]};
                    n_keys.push_back(new_key);
                    if(P.count(new_key) == 1){
                        if(P[new_key].prob >= G.thresh){
                            no_neighbors = false;
                        }
                    }      
                }
            }
        }

        All_Neighbors[current_key] = n_keys;  

        return no_neighbors;
    
    }else{

        for (int i = 0; i < All_Neighbors[current_key].size(); i++){
            vector<int> new_key = All_Neighbors[current_key][i];

            if(P.count(new_key) == 1){
                if(P[new_key].prob >= G.thresh){
                    return no_neighbors = false;
                }
            }  
        }
    }

    return no_neighbors;    
};

void HGBEES::Record_Data(string file_name, Grid G){
	ofstream myfile; myfile.open(file_name);
    for(int i = 0; i < n; i++){
        vector<int> current_key = keys[i]; 
        myfile << P[current_key].prob << " " << current_key[0] << " " << current_key[1] << " " << current_key[2] << endl;
    }  
    myfile.close();
};


#endif // HGBEES_H