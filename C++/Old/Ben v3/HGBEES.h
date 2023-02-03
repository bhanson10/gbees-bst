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

int state_count = 0; 
int key_count = 0; 

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

vector<int> UnshiftState(vector<int> state, int d){
    vector<int> unshift_state(d,0);
    for (int i = 0; i < d; i++){
        if(state[i] % 2 == 0){
            unshift_state[i] = state[i]/2;
        }else{
            unshift_state[i] = (state[i]+1)/-2;
        }
    }
    return unshift_state;
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

vector<int> CantorUnpair(int key, int d){
    vector<int> state(d,0);
    for (int i = 2; i <= d; i++){
        double z = key; double w = floor((sqrt(8*z+1)-1)/2); double t = (0.5)*(pow(w,2)+w); 
        double y = z - t; double x = w - y; state[d-i] = int(x); state[d-i+1] = int(y);
        key=x;
    }
    return state; 
};

int state_conversion(vector<int> state, Grid G){
    vector<int> shift_state = ShiftState(state, G.d);
    int key = CantorPair(shift_state);
    state_count++;
    return key;
};

vector<int> key_conversion(int key, Grid G){
    vector<int> shift_state = CantorUnpair(key,G.d);
    vector<int> state = UnshiftState(shift_state,G.d);
    key_count++;
    return state; 
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
    
    vector<int> state;

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                state = {i,j,k}; int key = state_conversion(state,G);
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[key].prob = exp(-4*x/2);
            }
        }
    }
    get_keys(); Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    int current_key; vector<int> current_state; vector<double> x(G.d,0); n = P.size(); 
    
    for (int l = 0; l < n; l++){
        current_key = keys[l]; 
        if((!P[current_key].active)&&(current_key!=-1)){
            current_state = key_conversion(current_key,G);
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
    get_keys();  n = P.size(); 

    for(int l = 0; l < n; l++){   // Check/Create Neighbors of Big Cells
        if (P[keys[l]].prob >= G.thresh){
            int current_key = keys[l]; vector<int> state = key_conversion(current_key, G);

            vector<double> current_v = P[current_key].v;

            vector<int> x_neighbors, y_neighbors, z_neighbors;

            if(current_v[0] > 0) x_neighbors = {state[0], state[0]+1}; else if (current_v[0] == 0) x_neighbors = {state[0]-1,state[0],state[0]+1}; else x_neighbors = {state[0]-1, state[0]};
            if(current_v[1] > 0) y_neighbors = {state[1], state[1]+1}; else if (current_v[1] == 0) y_neighbors = {state[1]-1,state[1],state[1]+1}; else y_neighbors = {state[1]-1, state[1]};
            if(current_v[2] > 0) z_neighbors = {state[2], state[2]+1}; else if (current_v[2] == 0) z_neighbors = {state[2]-1,state[2],state[2]+1}; else z_neighbors = {state[2]-1, state[2]};

            for (int i = 0; i < y_neighbors.size(); i++){
                for (int j = 0; j < x_neighbors.size(); j++){
                    for (int k = 0; k < z_neighbors.size(); k++){
                        vector<int> new_state = {x_neighbors[j], y_neighbors[i], z_neighbors[k]}; int new_key = state_conversion(new_state,G);
                        if(!P.count(new_key)){
                            P[new_key].prob = 0; P[new_key].active = 0; 
                        }                       
                    }
                }
            }               
        }
    }   
    n = P.size(); get_keys(); 

    for (int l = 0; l < n; l++){     // Remove small cells that do not neighbor big cells
        int current_key = keys[l];
        if(current_key!=-1){
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
        int l_key = keys[l];
        if(l_key != -1){
            vector<int> l_state = key_conversion(l_key, G);
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
                vector<int> l_state = key_conversion(l_key, G);
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
            vector<int> l_state = key_conversion(l_key, G);

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
    keys.clear();

    for (auto it = P.begin(); it != P.end(); it++) {
        keys.push_back(it->first);
    }
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, int current_key){
    bool neighbors = true; vector<int> state = key_conversion(current_key,G);

    for (int i = state[1]-1; i <= state[1]+1; i++){
        for (int j = state[0]-1; j <= state[0]+1; j++){
            for (int k = state[2]-1; k <= state[2]+1; k++){
                vector<int> new_state = {j, i, k}; int new_key = state_conversion(new_state,G);
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
    for(int i = 0; i < n; i++){
        int key = keys[i]; vector<int> state = key_conversion(key,G); 

        myfile << P[key].prob << " " << state[0] << " " << state[1] << " " << state[2] << endl;
    }  
    myfile.close();
};


#endif // HGBEES_H