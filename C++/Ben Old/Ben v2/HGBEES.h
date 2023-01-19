#ifndef HGBEES_H
#define HGBEES_H

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
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
    map<int, double> P;
    map<int, vector<double>> v;
    map<int, vector<double>> u;
    map<int, vector<double>> w;
    map<int, vector<double>> f;
    vector<int> keys;
    vector<double> values;
    int n;        

    void Initialize_D(Grid G,Lorenz3D Lor);
    void Initialize_vuw(Grid G,Lorenz3D Lor);
    void Modify_pointset(Grid G,Lorenz3D Lor);
    map<int, double> RHS_P(Grid G,Lorenz3D Lor);
    void get_keys_values();
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
    return key;
};

vector<int> key_conversion(int key, Grid G){
    vector<int> shift_state = CantorUnpair(key,G.d);
    vector<int> state = UnshiftState(shift_state,G.d);
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
    P[-1] = 0; v[-1] = zeros; u[-1] = zeros; w[-1] = zeros; f[-1] = zeros; 
    
    vector<int> state;

    for (int i = round((G.start[0]-2)/G.del); i <= round((G.start[0]+2)/G.del); i++){
        for (int j = round((G.start[1]-2)/G.del); j <= round((G.start[1]+2)/G.del); j++){
            for (int k = round((G.start[2]-2)/G.del); k <= round((G.start[2]+2)/G.del); k++){
                state = {i,j,k}; int key = state_conversion(state,G);
                double x = pow(i*G.del - G.start[0],2)+pow(j*G.del - G.start[1],2)+pow(k*G.del - G.start[2],2);
                P[key] = exp(-4*x/2);
            }
        }
    }
    get_keys_values(); Initialize_vuw(G,Lor); 
};

void HGBEES::Initialize_vuw(Grid G,Lorenz3D Lor){
    int current_key; vector<int> current_state; vector<double> x(G.d,0); n = P.size(); 
    
    for (int l = 0; l < n; l++){
        current_key = keys[l]; 
        if(current_key!=-1){
            current_state = key_conversion(current_key,G);
            for(int i = 0; i < G.d; i++){
                x[i] = G.del*current_state[i];
            }
            double v1 = Lor.sigma*(x[1]-(x[0]+G.xh));
            double v2 = -(x[1]+G.xh)-x[0]*x[2];
            double v3 = -Lor.b*(x[2]+G.xh)+x[0]*x[1]-Lor.b*Lor.r; 
            
            v[current_key] = {v1,v2,v3};
            u[current_key] = {min(v1,0.0),min(v2,0.0),min(v3,0.0)};
            w[current_key] = {max(v1,0.0),max(v2,0.0),max(v3,0.0)}; 
        }
    }
};

void HGBEES::Modify_pointset(Grid G,Lorenz3D Lor){
    get_keys_values();  n = P.size(); 

    for(int l = 0; l < n; l++){   // Check/Create Neighbors of Big Cells
        if (P[keys[l]] >= G.thresh){
            int current_key = keys[l]; vector<int> state = key_conversion(current_key, G);

            vector<double> current_v = v[current_key];

            vector<int> x_neighbors, y_neighbors, z_neighbors;

            if(current_v[0] > 0) x_neighbors = {state[0], state[0]+1}; else if (current_v[0] == 0) x_neighbors = {state[0]-1,state[0],state[0]+1}; else x_neighbors = {state[0]-1, state[0]};
            if(current_v[1] > 0) y_neighbors = {state[1], state[1]+1}; else if (current_v[1] == 0) y_neighbors = {state[1]-1,state[1],state[1]+1}; else y_neighbors = {state[1]-1, state[1]};
            if(current_v[2] > 0) z_neighbors = {state[2], state[2]+1}; else if (current_v[2] == 0) z_neighbors = {state[2]-1,state[2],state[2]+1}; else z_neighbors = {state[2]-1, state[2]};

            for (int i = 0; i < y_neighbors.size(); i++){
                for (int j = 0; j < x_neighbors.size(); j++){
                    for (int k = 0; k < z_neighbors.size(); k++){
                        vector<int> new_state = {x_neighbors[j], y_neighbors[i], z_neighbors[k]}; int new_key = state_conversion(new_state,G);
                        if(!P.count(new_key)){
                            P[new_key] = 0;
                        }                       
                    }
                }
            }               
        }
    }   
    n = P.size(); get_keys_values(); Initialize_vuw(G,Lor); 

    for (int l = 0; l < n; l++){     // Remove small cells that do not neighbor big cells
        int current_key = keys[l];
        if(current_key!=-1){
            if ((values[l] < G.thresh)&&(no_neighbors(G,Lor,current_key))){
                P.erase(current_key); v.erase(current_key); u.erase(current_key); w.erase(current_key); 
            }
        }
    }
    n = P.size(); get_keys_values(); double prob_sum;
    for (int i = 0; i < n; i++){
        P[keys[i]] = max(P[keys[i]],0.); prob_sum += P[keys[i]];
    }
    for (int i = 0; i < n; i++){
        P[keys[i]] = P[keys[i]]/prob_sum;
    }
    get_keys_values();
};

map<int, double> HGBEES::RHS_P(Grid G,Lorenz3D Lor){
    map<int, double> K; n = P.size(); get_keys_values(); vector<double> zeros(G.d,0);

    for(int l = 0; l < n; l++){ //Initializing K
        int l_key = keys[l];
        K[l_key] = 0;
    }   

    for(int l = 0; l < n; l++){ //Calculating Initial f
        int l_key = keys[l];
        if(l_key != -1){
            vector<int> l_state = key_conversion(l_key, G);
            vector<double> f_l = zeros; vector<double> u_l = u[l_key]; vector<double> w_l = w[l_key];    
            for(int q = 0; q < G.d; q++){
                vector<int> k_state = l_state; k_state[q] = k_state[q]+1; int k_key = state_conversion(k_state,G);
                if(P.count(k_key)==0) k_key = -1; 
                f_l[q] = w_l[q]*P[l_key]+u_l[q]*P[k_key];
            }
            f[l_key] = f_l;
        }
    }

    for(int q = 0; q < G.d; q++){
        for(int l = 0; l < n; l++){ // Calculating Total f
            int l_key = keys[l]; 
            if(l_key != -1){
                vector<int> l_state = key_conversion(l_key, G);
                vector<double> f_l = f[l_key]; vector<double> w_l = w[l_key];
                vector<int> i_state = l_state; i_state[q] = i_state[q]-1; 
                int i_key = state_conversion(i_state,G);
                if(P.count(i_key)==0) i_key = -1; 
                vector<double> f_i = f[i_key]; vector<double> w_i = w[i_key]; 
                vector<double> u_i = u[i_key]; vector<double> v_i = v[i_key];
                if ((P[l_key]>G.thresh)||(P[i_key]>G.thresh)){
                    double F = G.dt*(P[l_key]-P[i_key])/(2*G.del);
                    for(int e = 0; e < G.d; e++){
                        if (e!=q){
                            vector<int> j_state = l_state; j_state[e]=j_state[e]-1; 
                            int j_key = state_conversion(j_state,G); if(P.count(j_key)==0) j_key = -1; 
                            vector<double> f_j = f[j_key]; 
                            vector<double> w_j = w[j_key];
                            vector<double> u_j = u[j_key];
                            vector<int> p_state = i_state; p_state[e] = p_state[e]-1; 
                            int p_key = state_conversion(p_state,G); if(P.count(p_key)==0) p_key = -1;
                            vector<double> f_p = f[p_key]; 
                            vector<double> w_p = w[p_key];
                            vector<double> u_p = u[p_key];  
                            
                            f_l[e] = f_l[e]-w_l[e]*w_i[q]*F;
                            f_j[e] = f_j[e]-u_j[e]*w_i[q]*F; f[j_key] = f_j;
                            f_i[e] = f_i[e]-w_i[e]*u_i[q]*F; 
                            f_p[e] = f_p[e]-u_p[e]*u_i[q]*F; f[p_key] = f_p; 
                        }                       
                    }
                    f[l_key] = f_l; f[i_key] = f_i;

                    vector<int> i_i_state = i_state; i_i_state[q] = i_i_state[q]-1;
                    int i_i_key = state_conversion(i_i_state,G); if(P.count(i_i_key)==0) i_i_key = -1; 
                    vector<int> k_state = l_state; k_state[q] = k_state[q] + 1; 
                    int k_key = state_conversion(k_state,G); if(P.count(k_key)==0) k_key = -1; 
                    
                    double th,t; 
                    if (v_i[q]>0){
                        th = (P[i_key]-P[i_i_key])/(P[l_key]-P[i_key]);
                    }else{
                        th = (P[k_key]-P[l_key])/(P[l_key]-P[i_key]);
                    }

                    vector<double> f_i = f[i_key];

                    if(v_i[q]>=0){
                        t = v_i[q];
                    }else{
                        t = -v_i[q];
                    } 
                    
                    f_i[q] += t*(G.del/G.dt - t)*F*MC(th); f[i_key] = f_i;
                }
            }
        }
    }
    for(int l = 0; l < n; l++){
        int l_key = keys[l]; 
        if(l_key != -1){
            vector<int> l_state = key_conversion(l_key, G);
            vector<double> f_l = f[l_key]; 

            for(int q = 0; q < G.d; q++){
                vector<int> i_state = l_state; i_state[q] = i_state[q] - 1;
                int i_key = state_conversion(i_state,G); if(P.count(i_key)==0) i_key = -1;  
                vector<double> f_i = f[i_key]; 
                K[l_key] -= (f_l[q]-f_i[q])/G.del;  
            }
        }
    }
    
    return K;
};

void HGBEES::get_keys_values(){
    keys.clear(); values.clear();

    for (auto it = P.begin(); it != P.end(); it++) {
        keys.push_back(it->first);
        values.push_back(it->second);
    }
};

bool HGBEES::no_neighbors(Grid G, Lorenz3D Lor, int current_key){
    bool neighbors = true; vector<int> state = key_conversion(current_key,G);

    for (int i = state[1]-1; i <= state[1]+1; i++){
        for (int j = state[0]-1; j <= state[0]+1; j++){
            for (int k = state[2]-1; k <= state[2]+1; k++){
                vector<int> new_state = {j, i, k}; int new_key = state_conversion(new_state,G);
                if(P.count(new_key) == 1){
                    if(P[new_key] >= G.thresh){
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
        int key = keys[i]; vector<int> state = key_conversion(key,G); double value = values[i];

        myfile << value << " " << state[0] << " " << state[1] << " " << state[2] << endl;
    }  
    myfile.close();
};


#endif // HGBEES_H