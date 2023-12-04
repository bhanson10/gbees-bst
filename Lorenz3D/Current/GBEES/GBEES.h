#ifndef GBEES_H
#define GBEES_H

#include "BST.h"
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
class Traj{ // Properties of the trajectory
  public:            
    double sigma;        
    double b;
    double r;
    double T; 
};

class GBEES{       // GBEES Class
  public:
    BST P;
    TreeNode* dead; 
    int a_count = 0; 
    int tot_count = 1; 
    double cfl_min_dt = 1E10;

    void Initialize_D(Grid G,Traj Lor);
    void Initialize_vuw(Grid G, Traj Lor, TreeNode* r);
    void Initialize_ik_nodes(Grid G, TreeNode* r); 
    void Modify_pointset(Grid G, Traj Lor);
    TreeNode* create_neighbors(Grid G, TreeNode* r);
    TreeNode* prune_tree(Grid G, TreeNode* r);
    void mark_cells(Grid G, TreeNode* r);
    TreeNode* delete_cells(Grid G, TreeNode* r);
    void normalize_tree(Grid G, TreeNode* r);
    void get_sum(Grid G, TreeNode* r, double& prob_sum);
    void divide_sum(Grid G, TreeNode* r, double prob_sum);
    void RHS(Grid G, Traj Lor);
    void get_dcu(Grid G, TreeNode* r); 
    void update_ctu(Grid G, TreeNode* r); 
    void update_prob(Grid G, TreeNode* r); 
    void max_key(Grid G, TreeNode* r, uint64_t& key); 
    void check_CFL_condition(Grid G, TreeNode* r, double& cfl_min_dt); 
    void Record_Data(std::string file_name, Grid G, TreeNode* r, double t);
    void measurement_update(Measurement m, Grid G, TreeNode* r);
};
/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
uint64_t CantorPair(int state[], int n){
    uint64_t key; 
    if(n>2){
        int last = state[n-1]; int new_state[n-1];
        for (int i = 0; i < n - 1; i++){
            new_state[i] = state[i];
        }
        uint64_t x = CantorPair(new_state, n-1); int y = last;
        key = (0.5)*(x+y)*(x+y+1)+y;
    }else{
        int x = state[0]; int y = state[1];
        key = (0.5)*(x+y)*(x+y+1)+y;
    }
    return key; 
};

uint64_t RosenbergPair(int state[], int d, int m){
    uint64_t key; 
    if(d == 1){
        return key = state[0];
    }else{
        int new_state[d-1];
        for (int i = 0; i < d-1; i++){
            new_state[i] = state[i]; 
        }
        int new_m = *std::max_element(new_state, new_state + d-1);
        return key = RosenbergPair(new_state, d-1, new_m) + pow(m,d) + (m - state[d-1])*(pow(m+1, d-1) - pow(m,d-1)); 
    }
};

uint64_t state_conversion(std::array<int,DIM> state, Grid G){
    int shift_state[DIM];
    for (int i = 0; i < DIM; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    if(DIM == 1){
        uint64_t key = shift_state[0];
        return key;
    }else{
        uint64_t key = CantorPair(shift_state, DIM);
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
void GBEES::Initialize_D(Grid G, Traj Lor){
    Cell blank_c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .new_f = -1, .vuw_f = -1, .del_f = -1};
    TreeNode* dead_node = new TreeNode(-1, blank_c); dead_node->cell.i_nodes = {dead_node, dead_node, dead_node}; dead_node->cell.k_nodes = {dead_node, dead_node, dead_node};
    dead = dead_node; 

    std::array<int,DIM> current_state; uint64_t key; 
   for (int i = round((-3*G.std[0])/G.del[0]); i <= round((3*G.std[0])/G.del[0]); i++){
        for (int j = round((-3*G.std[1])/G.del[1]); j <= round((3*G.std[1])/G.del[1]); j++){
            for (int k = round((-3*G.std[2])/G.del[2]); k <= round((3*G.std[2])/G.del[2]); k++){
                current_state = {i,j,k}; key = state_conversion(current_state, G);

                double x = 0; 
                for(int q = 0; q < DIM; q++){
                    x += pow((current_state[q]*G.del[q]),2)/pow(G.std[q],2); 
                }
                
                Cell c = {.prob = exp(-x/2), .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = current_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                TreeNode* new_node = new TreeNode(key, c); P.root = P.insertRecursive(P.root, new_node);
            }
        }
    }

    Initialize_vuw(G,Lor,P.root);
    Initialize_ik_nodes(G,P.root); 
};

void GBEES::Initialize_vuw(Grid G, Traj Lor, TreeNode* r){
    if(r==NULL){ //Iterator over entire BST
        return;
    }
    Initialize_vuw(G,Lor,r->left);
    Initialize_vuw(G,Lor,r->right);  
    
    if(r->cell.new_f==0){
        std::array<double,DIM> x; 
        for(int i = 0; i < DIM; i++){
            x[i] = G.del[i]*r->cell.state[i]+G.epoch[i];
        }
        
        double v1 = Lor.sigma*(x[1]-(x[0]+G.xh[0]));
        double v2 = -(x[1]+G.xh[1])-x[0]*x[2];
        double v3 = -Lor.b*(x[2]+G.xh[2])+x[0]*x[1]-Lor.b*Lor.r; 
        r->cell.v = {v1,v2,v3};
        r->cell.u = {std::min(v1,0.0),std::min(v2,0.0),std::min(v3,0.0)};
        r->cell.w = {std::max(v1,0.0),std::max(v2,0.0),std::max(v3,0.0)}; 
        r->cell.new_f = 1;

        double sum = 0; 
        for(int q = 0; q < DIM; q++){
            sum += std::abs(r->cell.v[q])/G.del[q];
        }
        r->cell.cfl_dt = 1/sum; 
    }
};

void GBEES::Initialize_ik_nodes(Grid G,TreeNode* r){
    if (r==NULL){
        return; 
    }
    Initialize_ik_nodes(G,r->left);
    Initialize_ik_nodes(G,r->right);

    if(r->cell.vuw_f == 0){
        std::array<int, DIM> l_state = r->cell.state; 
        for(int q = 0; q < DIM; q++){
            //Initializing i, k nodes
            std::array<int,DIM> i_state = l_state; i_state[q] = i_state[q]-1; uint64_t i_key = state_conversion(i_state,G);
            std::array<int,DIM> k_state = l_state; k_state[q] = k_state[q]+1; uint64_t k_key = state_conversion(k_state,G);
            TreeNode* i_node = P.recursiveSearch(P.root, i_key); TreeNode* k_node = P.recursiveSearch(P.root, k_key); 
            
            if(i_node == NULL){
                i_node = dead; r->cell.i_nodes[q] = i_node; 
            }else{
                r->cell.i_nodes[q] = i_node; i_node->cell.k_nodes[q] = r; 
            }

            if(k_node == NULL){
                k_node = dead; r->cell.k_nodes[q] = k_node; 
            }else{
                r->cell.k_nodes[q] = k_node; k_node->cell.i_nodes[q] = r; 
            }      
        }
        r->cell.vuw_f = 1; 
    }
};

void GBEES::Modify_pointset(Grid G, Traj Lor){
    P.root = create_neighbors(G, P.root);
    Initialize_ik_nodes(G,P.root); 
    Initialize_vuw(G, Lor, P.root); 
    if(!P.isBalanced(P.root)){
        P.root = P.balance(P.root);
    }
};

/*
TreeNode* GBEES::create_neighbors(Grid G, TreeNode* r){
    if (r == NULL){
        return r;
    }
    create_neighbors(G, r->left);
    create_neighbors(G, r->right);

    if (r->cell.prob >= G.thresh){
        std::array<int,DIM> new_state = r->cell.state; uint64_t new_key; TreeNode* new_node; Cell c; 
        for (int q = 0; q < DIM; q++){

            // Checking Forward Faces
            if(r->cell.k_nodes[q] == dead){
                new_state[q] += 1; new_key = state_conversion(new_state,G);
                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                new_node = new TreeNode(new_key, c);
                P.root = P.insertRecursive(P.root, new_node); 
                // Checking Edges
                for (int e = 0; e < DIM; e++){
                    if(e != q){
                        new_state[e] += 1; new_key = state_conversion(new_state,G);
                        c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                        new_node = new TreeNode(new_key, c);
                        P.root = P.insertRecursive(P.root, new_node); 

                        new_state[e] -= 1; new_key = state_conversion(new_state,G);
                        c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                        new_node = new TreeNode(new_key, c);
                        P.root = P.insertRecursive(P.root, new_node); 
                    }
                }
            }else{
                new_state[q] += 1;
                for (int e = 0; e < DIM; e++){
                    if(e != q){

                        if(r->cell.k_nodes[q]->cell.k_nodes[e] == dead){
                            new_state[e] += 1; new_key = state_conversion(new_state,G);
                            c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                            new_node = new TreeNode(new_key, c);
                            P.root = P.insertRecursive(P.root, new_node); 
                        }
                    
                        if(r->cell.k_nodes[q]->cell.i_nodes[e] == dead){
                            new_state[e] -= 1; new_key = state_conversion(new_state,G);
                            c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                            new_node = new TreeNode(new_key, c);
                            P.root = P.insertRecursive(P.root, new_node); 
                        }
                    }
                }
            }

            // Checking Back Faces
            if(r->cell.i_nodes[q] == dead){
                new_state[q] -= 1; new_key = state_conversion(new_state,G);
                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                new_node = new TreeNode(new_key, c);
                P.root = P.insertRecursive(P.root, new_node); 

                // Checking Edges
                for (int e = 0; e < DIM; e++){
                    if(e != q){
                        new_state[e] += 1; new_key = state_conversion(new_state,G);
                        c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                        new_node = new TreeNode(new_key, c);
                        P.root = P.insertRecursive(P.root, new_node); 

                        new_state[e] -= 1; new_key = state_conversion(new_state,G);
                        c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                        new_node = new TreeNode(new_key, c);
                        P.root = P.insertRecursive(P.root, new_node); 
                    }
                }
            }else{
                new_state[q] -= 1;
                for (int e = 0; e < DIM; e++){
                    if(e != q){

                        if(r->cell.k_nodes[q]->cell.k_nodes[e] == dead){
                            new_state[e] += 1; new_key = state_conversion(new_state,G);
                            c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                            new_node = new TreeNode(new_key, c);
                            P.root = P.insertRecursive(P.root, new_node); 
                        }
                    
                        if(r->cell.k_nodes[q]->cell.i_nodes[e] == dead){
                            new_state[e] -= 1; new_key = state_conversion(new_state,G);
                            c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                            new_node = new TreeNode(new_key, c);
                            P.root = P.insertRecursive(P.root, new_node); 
                        }
                    }
                }
            }
        }
    }

    return P.root; 

};
*/

TreeNode* GBEES::create_neighbors(Grid G, TreeNode* r){
    if (r == NULL){
        return r;
    }
    create_neighbors(G, r->left);
    create_neighbors(G, r->right);

    //if ((r->cell.prob >= G.thresh)&&(!r->cell.center_f)){
    if (r->cell.prob >= G.thresh){
        std::array<double,DIM> current_v = r->cell.v; std::array<int,DIM> new_state = r->cell.state; uint64_t new_key; 
        for (int q = 0; q < DIM; q++){
            // Checking Faces
            if(current_v[q] > 0){
                if(r->cell.k_nodes[q] == dead){
                    new_state[q] += 1; new_key = state_conversion(new_state,G);
                    Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                    TreeNode* new_node = new TreeNode(new_key, c);
                    P.root = P.insertRecursive(P.root, new_node); 
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1; new_key = state_conversion(new_state,G);
                                Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                                TreeNode* new_node = new TreeNode(new_key, c);
                                P.root = P.insertRecursive(P.root, new_node); 
                            }else if (current_v[e] < 0){
                                new_state[e] -= 1; new_key = state_conversion(new_state,G);
                                Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                                TreeNode* new_node = new TreeNode(new_key, c);
                                P.root = P.insertRecursive(P.root, new_node); 
                            }
                        }
                    }
                }
            // Checking Faces
            }else if(current_v[q] < 0){
                if(r->cell.i_nodes[q] == dead){
                    new_state[q] -= 1; new_key = state_conversion(new_state,G);
                    Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                    TreeNode* new_node = new TreeNode(new_key, c);
                    P.root = P.insertRecursive(P.root, new_node); 

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1; new_key = state_conversion(new_state,G);
                                Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                                TreeNode* new_node = new TreeNode(new_key, c);
                                P.root = P.insertRecursive(P.root, new_node); 
                            }else if (current_v[e] < 0){
                                new_state[e] -= 1; new_key = state_conversion(new_state,G);
                                Cell c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .state = new_state, .new_f = 0, .vuw_f = 0, .del_f = 0};
                                TreeNode* new_node = new TreeNode(new_key, c);
                                P.root = P.insertRecursive(P.root, new_node); 
                            }
                        }
                    }
                }
            }
        }
    }

    return P.root; 

};

TreeNode* GBEES::prune_tree(Grid G, TreeNode* r){
    mark_cells(G, r); 
    P.root = delete_cells(G, r); 

    return P.root; 
};

void GBEES::mark_cells(Grid G, TreeNode* r){
    if (r==NULL){
        return;
    }
    mark_cells(G, r->left);
    mark_cells(G, r->right);

    r->cell.vuw_f = 0; bool DELETE = true; 
    if (r->cell.prob < G.thresh){
        
        for(int q = 0; q < DIM; q++){
            // Looking at Backwards Node

            if(r->cell.i_nodes[q]!=dead){
                //if (r->cell.i_nodes[q]->cell.prob >= G.thresh){
                if ((r->cell.i_nodes[q]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.prob >= G.thresh)){
                    DELETE = false; 
                    break; 
                }else{
                    for (int e = 0; e < DIM; e++){
                        if(e!=q){
                            if ((r->cell.i_nodes[q]->cell.i_nodes[e]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.v[e]>0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh)){
                            //if (r->cell.i_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh){
                                DELETE = false; 
                                break; 
                            }

                            if ((r->cell.i_nodes[q]->cell.k_nodes[e]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.k_nodes[e]->cell.v[e]<0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh)){
                            //if (r->cell.i_nodes[q]->cell.k_nodes[e]->cell.prob >= G.thresh){
                                DELETE = false; 
                                break; 
                            }
                        }
                    }
                }
            }
            // Looking at Forwards Node
            if(r->cell.k_nodes[q]!=dead){
                //if (r->cell.k_nodes[q]->cell.prob >= G.thresh){
                if ((r->cell.k_nodes[q]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.prob >= G.thresh)){
                    DELETE = false; 
                    break; 
                }else{
                    for (int e = 0; e < DIM; e++){
                        if(e!=q){
                            //if (r->cell.k_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh){
                            if ((r->cell.k_nodes[q]->cell.i_nodes[e]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.i_nodes[e]->cell.v[e]>0)&&(r->cell.k_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh)){
                                DELETE = false; 
                                break; 
                            }

                            //if (r->cell.k_nodes[q]->cell.k_nodes[e]->cell.prob >= G.thresh){
                            if ((r->cell.k_nodes[q]->cell.k_nodes[e]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.k_nodes[e]->cell.v[e]<0)&&(r->cell.k_nodes[q]->cell.k_nodes[e]->cell.prob >= G.thresh)){
                                DELETE = false; 
                                break; 
                            }
                        }
                    }
                }
            }
        }

        if(DELETE){
            r->cell.del_f = 1; 
        }
    }
};

/*
void GBEES::mark_cells(Grid G, TreeNode* r){
    if (r==NULL){
        return;
    }
    mark_cells(G, r->left);
    mark_cells(G, r->right);

    r->cell.vuw_f = 0; bool DELETE = true; 
    if (r->cell.prob < G.thresh){
        
        for(int q = 0; q < DIM; q++){
            // Looking at Backwards Node
            if(r->cell.i_nodes[q]!=dead){
                if ((r->cell.i_nodes[q]->cell.v[q] > 0)&&(r->cell.i_nodes[q]->cell.prob >= G.thresh)){
                    DELETE = false; 
                    break; 
                }
            }
            // Looking at Forwards Node
            if(r->cell.k_nodes[q]!=dead){
                if ((r->cell.v[q] < 0)&&(r->cell.k_nodes[q]->cell.prob >= G.thresh)){
                    DELETE = false; 
                    break; 
                }
            }
        }

        if(DELETE){
            r->cell.del_f = 1; 
        }
    }
};
*/

TreeNode* GBEES::delete_cells(Grid G, TreeNode* r){
    if (r==NULL){
        return r;
    }
    delete_cells(G, r->left);
    delete_cells(G, r->right);

    if (r->cell.del_f == 1){
        P.root = P.deleteNode(P.root, r->key);
    }

    return P.root; 
};

void GBEES::normalize_tree(Grid G, TreeNode* r){
    double prob_sum = 0; a_count = 0; tot_count = 1; 

    get_sum(G, r, prob_sum);
    divide_sum(G,r,prob_sum);
};

void GBEES::get_sum(Grid G, TreeNode* r, double& prob_sum){
    if (r==NULL){
        return;
    }
    get_sum(G, r->left, prob_sum);
    get_sum(G, r->right, prob_sum);

    if(r->cell.prob >= G.thresh){
        a_count++; 
    }
    prob_sum += r->cell.prob; 
    tot_count++; 
};
void GBEES::divide_sum(Grid G, TreeNode* r, double prob_sum){
    if (r==NULL){
        return;
    }
    divide_sum(G, r->left, prob_sum);
    divide_sum(G, r->right, prob_sum);

    r->cell.prob /= prob_sum; 
};

void GBEES::RHS(Grid G, Traj Lor){
    get_dcu(G, P.root); 
    update_ctu(G, P.root); 
};

void GBEES::get_dcu(Grid G, TreeNode* r){
    if (r==NULL){
        return; 
    }
    get_dcu(G,r->left);
    get_dcu(G,r->right);

    if(r->key != -1){
        r->cell.dcu = 0; r->cell.ctu = {0}; std::array <double, DIM> mu = {0, 0, 0};
        for(int q = 0; q < DIM; q++){
            TreeNode* i_node = r->cell.i_nodes[q]; 
            
            double dcu_p = (r->cell.w[q] * r->cell.prob + r->cell.u[q] * r->cell.k_nodes[q]->cell.prob) - (mu[q])*(r->cell.k_nodes[q]->cell.prob-r->cell.prob)/G.del[q];
            double dcu_m = (i_node->cell.w[q]*i_node->cell.prob + i_node->cell.u[q]*r->cell.prob)  - (mu[q])*(r->cell.prob-i_node->cell.prob)/G.del[q];
            r->cell.dcu -= (G.dt/G.del[q])*(dcu_p-dcu_m); 
        }
    }
};

void GBEES::update_ctu(Grid G, TreeNode* r){
    if (r==NULL){
        return;
    }
    update_ctu(G,r->left);
    update_ctu(G,r->right);

    if(r->key != -1){ 
        for(int a = 0; a < DIM; a++){
            TreeNode* i_node = r->cell.i_nodes[a]; 
            TreeNode* j_node; TreeNode* p_node;
            if(i_node->key!=-1){
                double F = G.dt*(r->cell.prob-i_node->cell.prob)/(2*G.del[a]); 
                for(int b = 0; b < DIM; b++){
                    if (b!=a){
                        j_node = r->cell.i_nodes[b]; 
                        p_node = i_node->cell.i_nodes[b]; 

                        r->cell.ctu[b]      -= i_node->cell.w[a] * r->cell.w[b] * F;
                        j_node->cell.ctu[b] -= i_node->cell.w[a] * j_node->cell.u[b] * F;
                        i_node->cell.ctu[b] -= i_node->cell.u[a] * i_node->cell.w[b] * F;
                        p_node->cell.ctu[b] -= i_node->cell.u[a] * p_node->cell.u[b] * F; 
                    }                 
                }
                
                // High-Resolution Correction Terms
                double th; 
                if (i_node->cell.v[a]>0){
                    TreeNode* i_i_node = i_node->cell.i_nodes[a];
                    th = (i_node->cell.prob-i_i_node->cell.prob)/(r->cell.prob-i_node->cell.prob);
                }else{
                    th = (r->cell.k_nodes[a]->cell.prob-r->cell.prob)/(r->cell.prob-i_node->cell.prob);
                }
                
                i_node->cell.ctu[a] += std::abs(i_node->cell.v[a])*(G.del[a]/G.dt - std::abs(i_node->cell.v[a]))*F*MC(th); 
            }
        }
    }
};

void GBEES::update_prob(Grid G, TreeNode* r){
    if (r == NULL){
        return; 
    }

    update_prob(G, r->left);
    update_prob(G, r->right);
    
    if(r->key != -1){
        r->cell.prob += r->cell.dcu; 
        for(int q = 0; q < DIM; q++){
            r->cell.prob -= (G.dt/G.del[q])*(r->cell.ctu[q]-r->cell.i_nodes[q]->cell.ctu[q]); 
        }
    }
};

void GBEES::Record_Data(std::string file_name, Grid G, TreeNode* r, double t){
	std::ofstream myfile; myfile.open(file_name);
    myfile << t << std::endl; 
    P.writeFile(myfile, G, P.root);
    myfile.close(); 
};

void GBEES::max_key(Grid G, TreeNode* r, uint64_t& key){
    if (r == NULL){
        return; 
    }

    max_key(G, r->left, key);
    max_key(G, r->right, key);

    key = std::max(key, r->key); 
};

void GBEES::measurement_update(Measurement m, Grid G, TreeNode* r){

    if (r == NULL){
        return; 
    }

    measurement_update(m, G, r->left);
    measurement_update(m, G, r->right);

    double x = 0; 
    for(int q = 0; q < DIM; q++){
        if (m.std[q] != 0){
            x += pow((r->cell.state[q]*G.del[q]-m.mean[q]),2)/pow(m.std[q],2); 
        }
    }

    double prob = exp(-x/2); r->cell.prob *= prob; 
};

void GBEES::check_CFL_condition(Grid G, TreeNode* r, double& cfl_min_dt){
    if (r == NULL){
        return; 
    }

    check_CFL_condition(G, r->left, cfl_min_dt);
    check_CFL_condition(G, r->right, cfl_min_dt);

    if(r->key != -1){
        cfl_min_dt = std::min(cfl_min_dt,r->cell.cfl_dt);
    }
};

#endif // GBEES_H