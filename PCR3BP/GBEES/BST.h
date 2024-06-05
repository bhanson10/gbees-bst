//
// Created by Benjamin Hanson on 5/13/24.
//

#ifndef GBEES_BST_H
#define GBEES_BST_H

#include <iostream>
#include <cwchar>
#include <array>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>
#include <sstream>
#include <armadillo>
const int DIM = 4;
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
class TreeNode;

class Cell{  // Cell - data connected to a given cell
public:
    double prob;
    std::array <double, DIM> v;
    std::array <double, DIM> u;
    std::array <double, DIM> w;
    std::array <double, DIM> ctu;
    std::array <int, DIM> state;
    std::array<TreeNode*, DIM> i_nodes;
    std::array<TreeNode*, DIM> k_nodes;
    double dcu;
    double cfl_dt;
    int new_f;
    int ik_f;
    int del_f;
};

class Grid{ // Properties of the grid
public:
    bool DIFF_B{};
    double thresh{};
    double dt{};
    std::array<double,DIM> epoch{};
    std::array<double, DIM> del{};
    arma::mat cov;
    std::array <double, DIM> diff{};

    explicit Grid(int n) : cov(n, n, arma::fill::zeros) {}
};

class Traj{
public:
    double mu;
    double T;
};

class TreeNode{ // TreeNode - Node in a Binary search tree
public:
    uint64_t key;
    Cell cell{};
    TreeNode* left;
    TreeNode* right;
    TreeNode(uint64_t k, Cell c);
};

class Measurement{
public:
    std::array<double, DIM> mean;
    std::array<double, DIM> std;
};

class BST{ // Binary Search Tree - data structure storing information
public:
    TreeNode* dead;
    TreeNode* root;
    int a_count = 0;
    int tot_count = 1;
    double cfl_min_dt = 1E10;
    BST();
    TreeNode* insert_recursive(TreeNode* r, TreeNode* new_node);
    TreeNode* recursive_search(TreeNode* r, uint64_t k);
    static TreeNode* min_value_node(TreeNode* node);
    TreeNode* delete_node(TreeNode* r, uint64_t k);
    int get_height(TreeNode* r);
    int get_difference(TreeNode* r);
    static TreeNode* rr_rotate(TreeNode* parent);
    static TreeNode* ll_rotate(TreeNode* parent);
    static TreeNode* lr_rotate(TreeNode* parent);
    static TreeNode* rl_rotate(TreeNode* parent);
    TreeNode* balance(TreeNode* r);
    bool is_balanced(TreeNode* r);
    void write_file(std::ofstream& myfile, const Grid& G, TreeNode* r);
    void initialize_grid(const Grid& G, Traj lyap, Measurement m);
    void initialize_vuw(const Grid& G, Traj lyap, TreeNode* r);
    void initialize_ik_nodes(TreeNode* r);
    void initialize_ik_nodes_single(TreeNode* r);
    void grow_tree_1(const Grid& G, Traj lyap);
    void create_neighbors_1(const Grid& G, TreeNode* r);
    void grow_tree_2(const Grid& G, Traj lyap);
    void create_neighbors_2(const Grid& G, TreeNode* r);
    void prune_tree(const Grid& G, TreeNode* r);
    void mark_cells(const Grid& G, TreeNode* r);
    void delete_cells(const Grid& G, TreeNode* r);
    void normalize_tree(const Grid& G, TreeNode* r);
    void get_sum(const Grid& G, TreeNode* r, double& prob_sum);
    void divide_sum(const Grid& G, TreeNode* r, double prob_sum);
    void godunov_method(const Grid& G);
    void get_dcu(const Grid& G, TreeNode* r);
    void update_ctu(const Grid& G, TreeNode* r);
    void update_prob(const Grid& G, TreeNode* r);
    void max_key(const Grid& G, TreeNode* r, uint64_t& key);
    void check_cfl_condition(const Grid& G, TreeNode* r);
    void record_data(const std::string& file_name, const Grid& G, double t);
    void measurement_update(Measurement m, const Grid& G, TreeNode* r);
};
/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
uint64_t rosenberg_pair(const int* state, int d, int m){
    if(d == 1){
        return state[0];
    }

    int new_state[d-1];
    for (int i = 0; i < d-1; i++){
        new_state[i] = state[i];
    }
    int new_m = *std::max_element(new_state, new_state + d-1);
    return rosenberg_pair(new_state, d-1, new_m) + uint64_t(pow(m,d)) + (m - state[d-1])*uint64_t((pow(m+1, d-1)) - pow(m,d-1));
}

uint64_t state_conversion(const int* state){
    int* shift_state = new int[DIM]; int m; uint64_t key;
    for (int i = 0; i < DIM; i++){
        if(state[i]<0){
            shift_state[i] = -2*state[i]-1;
        }else{
            shift_state[i] = 2*state[i];
        }
    }
    m = *std::max_element(shift_state, shift_state + DIM);
    key = rosenberg_pair(shift_state, DIM, m);
    return key;
}

double mc(double th){
    double phi;
    phi = std::max(0.0,std::min({(1+th)/2,2.0,2*th}));
    return phi;
}

int get_size(TreeNode* r){
    if (r == nullptr){
        return 0;
    }

    return 1 + get_size(r->left) + get_size(r->right);
}
/*==============================================================================
Member Function DEFINITIONS
==============================================================================*/
TreeNode::TreeNode(uint64_t k, Cell c){ // Initialization of Node
    key = k;
    cell = c;
    left = nullptr; right = nullptr;
}

BST::BST(){ // Initialization of BST
    root = nullptr;
    dead = nullptr;
}

TreeNode* BST::insert_recursive(TreeNode* r, TreeNode* new_node){ // Inserting a node into BST
    if (r==nullptr){
        r=new_node;
        return r;
    }
    if (new_node->key < r->key){
        r->left = insert_recursive(r->left, new_node);
    }else if (new_node->key > r->key){
        r->right = insert_recursive(r->right, new_node);
    }else{
        return r;
    }
    return r;
}

TreeNode* BST::recursive_search(TreeNode* r, uint64_t k){ // Searching to see if a key exists in BST
    if ((r==nullptr)||(r->key == k)){
        return r;
    }else if (k < r->key){
        return recursive_search(r->left, k);
    }else{
        return recursive_search(r->right,k);
    }
}

TreeNode* BST::min_value_node(TreeNode* node){ // Getting the minimum value of a node for deletion
    TreeNode* current = node;
    while(current->left != nullptr){
        current = current->left;
    }
    return current;
}

TreeNode* BST::delete_node(TreeNode* r, uint64_t k){ // Deleting a node from the BST
    if(r==nullptr){
        return nullptr;
    }else if(k < r->key){
        r->left = delete_node(r->left, k);
    }else if(k > r->key){
        r->right = delete_node(r->right,k);
    }else{
        if(r->left == nullptr){
            TreeNode* temp = r->right;
            delete r;
            return temp;
        } else if (r->right == nullptr){
            TreeNode* temp = r->left;
            delete r;
            return temp;
        }else{
            TreeNode* temp = min_value_node(r->right);
            r->key = temp->key; r->cell = temp->cell;
            r->right = delete_node(r->right, temp->key);
        }
    }
    return r;
}

int BST::get_height(TreeNode* r){
    if (r == nullptr){
        return 0;
    }

    return 1 + std::max(get_height(r->left), get_height(r->right));
}

int BST::get_difference(TreeNode* r){
    int l_height = get_height(r->left);
    int r_height = get_height(r->right);
    int b_factor = l_height - r_height;
    return b_factor;
}

TreeNode* BST::rr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = t->left;
    t->left = parent;
    return t;
}

TreeNode* BST::ll_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = t->right;
    t->right = parent;
    return t;
}

TreeNode* BST::lr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = rr_rotate(t);
    return ll_rotate(parent);
}

TreeNode* BST::rl_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = ll_rotate(t);
    return rr_rotate(parent);
}

TreeNode* BST::balance(TreeNode* r){
    int bal_factor = get_difference(r);
    while(abs(bal_factor) > 1){
        if (bal_factor > 1){
            if (get_difference(r->left) > 0){
                r = ll_rotate(r);
            }else{
                r = lr_rotate(r);
            }
        }else if (bal_factor < -1){
            if (get_difference(r->right) > 0){
                r = rl_rotate(r);
            }else{
                r = rr_rotate(r);
            }
        }
        bal_factor = get_difference(r);
    }
    return r;
}

bool BST::is_balanced(TreeNode* r){
    int bal_factor = get_difference(r);
    if(abs(bal_factor) > 1){
        return false;
    }else{
        return true;
    }
}

void BST::write_file(std::ofstream& myfile, const Grid& G, TreeNode* r){
    if (r == nullptr){
        return;
    }

    write_file(myfile, G, r->left);
    write_file(myfile, G, r->right);

    if (r->cell.prob >= G.thresh){
        myfile << r->cell.prob;
        for(int i = 0; i < DIM; i++){
            myfile << " " << G.del[i]*r->cell.state[i]+G.epoch[i];
        }
        myfile << std::endl;
    }
}

void BST::initialize_grid(const Grid& G, Traj lyap, Measurement m){
    Cell blank_c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .dcu = 0, .new_f = -1, .ik_f = -1, .del_f = -1};
    auto* dead_node = new TreeNode(-1, blank_c); dead_node->cell.i_nodes = {dead_node, dead_node, dead_node, dead_node}; dead_node->cell.k_nodes = {dead_node, dead_node, dead_node, dead_node};
    dead = dead_node;

    std::array<int,DIM> current_state{}; uint64_t key;
    for (int i = int(round((-3*m.std[0])/G.del[0])); i <= int(round((3*m.std[0])/G.del[0])); i++){
        for (int j = int(round((-3*m.std[1])/G.del[1])); j <= int(round((3*m.std[1])/G.del[1])); j++){
            for (int k = int(round((-3*m.std[2])/G.del[2])); k <= int(round((3*m.std[2])/G.del[2])); k++){
                for (int l = int(round((-3*m.std[3])/G.del[3])); l <= int(round((3*m.std[3])/G.del[3])); l++){

                    current_state = {i,j,k,l}; key = state_conversion(current_state.data());
                    arma::vec current_state_vec = {i*G.del[0], j*G.del[1], k*G.del[2], l*G.del[3]};

                    double x = as_scalar(arma::trans(current_state_vec)*arma::inv(G.cov)*(current_state_vec));

                    Cell c = {.prob = exp(-x/2), .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = current_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};

                    if(c.prob >= G.thresh){
                        a_count++;
                    }
                    tot_count++;
                    auto* new_node = new TreeNode(key, c); root = insert_recursive(root, new_node);
                }
            }
        }
    }

    initialize_vuw(G,lyap,root);
    initialize_ik_nodes(root);
}

void BST::initialize_vuw(const Grid& G, Traj lyap, TreeNode* r){
    if(r==nullptr){ //Iterator over entire BST
        return;
    }
    initialize_vuw(G,lyap,r->left);
    initialize_vuw(G,lyap,r->right);

    if(r->cell.new_f==0){
        std::array<double,DIM> x{};
        for(int i = 0; i < DIM; i++){
            x[i] = G.del[i]*r->cell.state[i]+G.epoch[i];
        }

        double r1 = pow(pow(x[0]+lyap.mu,2)+pow(x[1],2),1.5);
        double r2 = pow(pow(x[0]-1+lyap.mu,2)+pow(x[1],2),1.5);

        double v1 = x[2];
        double v2 = x[3];
        double v3 = 2*x[3]+x[0]-(lyap.mu*(x[0]-1+lyap.mu)/r2)-((1-lyap.mu)*(x[0]+lyap.mu)/r1);
        double v4 = -2*x[2]+x[1]-(lyap.mu*x[1]/r2)-((1-lyap.mu)*x[1]/r1);

        r->cell.v = {v1,v2,v3,v4};
        r->cell.u = {std::min(v1,0.0),std::min(v2,0.0),std::min(v3,0.0),std::min(v4,0.0)};
        r->cell.w = {std::max(v1,0.0),std::max(v2,0.0),std::max(v3,0.0),std::max(v4,0.0)};
        r->cell.new_f = 1;

        double sum = 0;
        for(int q =0; q < DIM; q++){
            sum += std::abs(r->cell.v[q])/G.del[q];
        }

        double C_max;
        if(G.DIFF_B){
            C_max = 0.5;
        }else{
            C_max = 1;
        }

        r->cell.cfl_dt = C_max/sum;
    }
}

void BST::initialize_ik_nodes(TreeNode* r){
    if (r==nullptr){
        return;
    }
    initialize_ik_nodes(r->left);
    initialize_ik_nodes(r->right);

    initialize_ik_nodes_single(r);
}

void BST::initialize_ik_nodes_single(TreeNode* r){
    if(r->cell.ik_f == 0){
        std::array<int, DIM> l_state = r->cell.state;
        for(int q =0; q < DIM; q++){
            //Initializing i, k nodes
            std::array<int,DIM> i_state = l_state; i_state[q] = i_state[q]-1; uint64_t i_key = state_conversion(i_state.data());
            std::array<int,DIM> k_state = l_state; k_state[q] = k_state[q]+1; uint64_t k_key = state_conversion(k_state.data());
            TreeNode* i_node = recursive_search(root, i_key); TreeNode* k_node = recursive_search(root, k_key);

            if(i_node == nullptr){
                i_node = dead; r->cell.i_nodes[q] = i_node;
            }else{
                r->cell.i_nodes[q] = i_node; i_node->cell.k_nodes[q] = r;
            }

            if(k_node == nullptr){
                k_node = dead; r->cell.k_nodes[q] = k_node;
            }else{
                r->cell.k_nodes[q] = k_node; k_node->cell.i_nodes[q] = r;
            }
        }
        r->cell.ik_f = 1;
    }
}

void BST::grow_tree_1(const Grid& G, Traj lyap){
    create_neighbors_1(G, root);
    initialize_ik_nodes(root);
    initialize_vuw(G, lyap, root);
    if(!is_balanced(root)){
        root = balance(root);
    }
}

void BST::grow_tree_2(const Grid& G, Traj lyap){
    create_neighbors_2(G, root);
    initialize_vuw(G, lyap, root);
    if(!is_balanced(root)){
        root = balance(root);
    }
}

void BST::create_neighbors_1(const Grid& G, TreeNode* r){
    if (r == nullptr){
        return;
    }
    create_neighbors_1(G, r->left);
    create_neighbors_1(G, r->right);

    if (r->cell.prob >= G.thresh){
        Cell c{};
        TreeNode* new_node;
        std::array<double,DIM> current_v = r->cell.v;
        std::array<int,DIM> current_state = r->cell.state; std::array<int,DIM> new_state{}; uint64_t new_key;
        for(int q =0; q < DIM; q++){
            new_state = current_state;
            // Checking Forward Faces
            if(current_v[q] > 0){
                if(r->cell.k_nodes[q] == dead){
                    new_state[q] += 1;
                    new_key = state_conversion(new_state.data());
                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                    new_node = new TreeNode(new_key, c);
                    root = insert_recursive(root, new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                root = insert_recursive(root, new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                std::array<int, DIM> l_state = new_node->cell.state;
                                root = insert_recursive(root, new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                if(r->cell.k_nodes[q]->cell.k_nodes[e] == dead){
                                    new_state[e] += 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.k_nodes[q]->cell.i_nodes[e] == dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);

                                }
                            }
                        }
                    }
                }
                // Checking Backward Faces
            }else if(current_v[q] < 0){
                if(r->cell.i_nodes[q] == dead){
                    new_state[q] -= 1;
                    new_key = state_conversion(new_state.data());
                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                    new_node = new TreeNode(new_key, c);
                    root = insert_recursive(root, new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0,  .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                root = insert_recursive(root, new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                root = insert_recursive(root, new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                if(r->cell.i_nodes[q]->cell.k_nodes[e] == dead){
                                    new_state[e] += 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.i_nodes[q]->cell.i_nodes[e] == dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void BST::create_neighbors_2(const Grid& G, TreeNode* r){
    if (r == nullptr){
        return;
    }
    create_neighbors_2(G, r->left);
    create_neighbors_2(G, r->right);

    if (r->cell.prob >= G.thresh){
        Cell c{};
        TreeNode* new_node;
        std::array<double,DIM> current_v = r->cell.v;
        std::array<int,DIM> current_state = r->cell.state; std::array<int,DIM> new_state{}; uint64_t new_key;
        for(int q =0; q < DIM; q++){
            new_state = current_state;
            // Checking Forward Faces
            if(current_v[q] > 0){
                if(r->cell.k_nodes[q] == dead){
                    new_state[q] += 1;
                    new_key = state_conversion(new_state.data());
                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                    new_node = new TreeNode(new_key, c);
                    root = insert_recursive(root, new_node);
                    initialize_ik_nodes_single(new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                root = insert_recursive(root, new_node);
                                initialize_ik_nodes_single(new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                std::array<int, DIM> l_state = new_node->cell.state;
                                root = insert_recursive(root, new_node);
                                initialize_ik_nodes_single(new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                if(r->cell.k_nodes[q]->cell.k_nodes[e] == dead){
                                    new_state[e] += 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);
                                    initialize_ik_nodes_single(new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.k_nodes[q]->cell.i_nodes[e] == dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);
                                    initialize_ik_nodes_single(new_node);

                                }
                            }
                        }
                    }
                }
                // Checking Backward Faces
            }else if(current_v[q] < 0){
                if(r->cell.i_nodes[q] == dead){
                    new_state[q] -= 1;
                    new_key = state_conversion(new_state.data());
                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                    new_node = new TreeNode(new_key, c);
                    root = insert_recursive(root, new_node);
                    initialize_ik_nodes_single(new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0,  .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                root = insert_recursive(root, new_node);
                                initialize_ik_nodes_single(new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state.data());
                                c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                new_node = new TreeNode(new_key, c);
                                root = insert_recursive(root, new_node);
                                initialize_ik_nodes_single(new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        new_state = current_state; new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                if(r->cell.i_nodes[q]->cell.k_nodes[e] == dead){
                                    new_state[e] += 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);
                                    initialize_ik_nodes_single(new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.i_nodes[q]->cell.i_nodes[e] == dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state.data());
                                    c = {.prob = 0, .v = {0}, .u = {0}, .w = {0}, .ctu = {0}, .state = new_state, .dcu = 0, .new_f = 0, .ik_f = 0, .del_f = 0};
                                    new_node = new TreeNode(new_key, c);
                                    root = insert_recursive(root, new_node);
                                    initialize_ik_nodes_single(new_node);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void BST::prune_tree(const Grid& G, TreeNode* r){
    mark_cells(G, r);
    delete_cells(G, r);
    normalize_tree(G, root);
    initialize_ik_nodes(root);
}

void BST::mark_cells(const Grid& G, TreeNode* r){
    if (r==nullptr){
        return;
    }
    mark_cells(G, r->left);
    mark_cells(G, r->right);

    r->cell.ik_f = 0; bool DELETE = true;
    if (r->cell.prob < G.thresh){

        for(int q =0; q < DIM; q++){
            // Looking at Backwards Node

            if(r->cell.i_nodes[q]!=dead){
                if ((r->cell.i_nodes[q]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.prob >= G.thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int e = 0; e < DIM; e++){
                        if(e!=q){
                            if ((r->cell.i_nodes[q]->cell.i_nodes[e]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.v[e]>0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh)){
                                DELETE = false;
                                break;
                            }

                            if ((r->cell.i_nodes[q]->cell.k_nodes[e]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.k_nodes[e]->cell.v[e]<0)&&(r->cell.i_nodes[q]->cell.k_nodes[e]->cell.prob >= G.thresh)){
                                DELETE = false;
                                break;
                            }
                        }
                    }
                }
            }
            // Looking at Forwards Node
            if(r->cell.k_nodes[q]!=dead){
                if ((r->cell.k_nodes[q]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.prob >= G.thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int e = 0; e < DIM; e++){
                        if(e!=q){
                            if ((r->cell.k_nodes[q]->cell.i_nodes[e]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.i_nodes[e]->cell.v[e]>0)&&(r->cell.k_nodes[q]->cell.i_nodes[e]->cell.prob >= G.thresh)){
                                DELETE = false;
                                break;
                            }

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
}

void BST::delete_cells(const Grid& G, TreeNode* r){
    if (r==nullptr){
        return;
    }
    delete_cells(G, r->left);
    delete_cells(G, r->right);

    if (r->cell.del_f == 1){
        root = delete_node(root, r->key);
    }
}

void BST::normalize_tree(const Grid& G, TreeNode* r){
    double prob_sum = 0; a_count = 0; tot_count = 1;

    get_sum(G, r, prob_sum);
    divide_sum(G,r,prob_sum);
}

void BST::get_sum(const Grid& G, TreeNode* r, double& prob_sum){
    if (r==nullptr){
        return;
    }
    get_sum(G, r->left, prob_sum);
    get_sum(G, r->right, prob_sum);

    if(r->cell.prob >= G.thresh){
        a_count++;
    }
    prob_sum += r->cell.prob;
    tot_count++;
}

void BST::divide_sum(const Grid& G, TreeNode* r, double prob_sum){
    if (r==nullptr){
        return;
    }
    divide_sum(G, r->left, prob_sum);
    divide_sum(G, r->right, prob_sum);

    r->cell.prob /= prob_sum;
}

void BST::godunov_method(const Grid& G){
    get_dcu(G, root);
    update_ctu(G, root);
}

void BST::get_dcu(const Grid& G, TreeNode* r){
    if (r==nullptr){
        return;
    }
    get_dcu(G,r->left);
    get_dcu(G,r->right);

    r->cell.dcu = 0; r->cell.ctu = {0};
    for(int q =0; q < DIM; q++){
        TreeNode* i_node = r->cell.i_nodes[q];

        double dcu_p, dcu_m;
        if(G.DIFF_B){
            dcu_p = (r->cell.w[q] * r->cell.prob + r->cell.u[q] * r->cell.k_nodes[q]->cell.prob) - (G.diff[q])*(r->cell.k_nodes[q]->cell.prob-r->cell.prob)/G.del[q];
            dcu_m = (i_node->cell.w[q]*i_node->cell.prob + i_node->cell.u[q]*r->cell.prob)  - (G.diff[q])*(r->cell.prob-r->cell.i_nodes[q]->cell.prob)/G.del[q];
        }else{
            dcu_p = (r->cell.w[q] * r->cell.prob + r->cell.u[q] * r->cell.k_nodes[q]->cell.prob);
            dcu_m = (i_node->cell.w[q]*i_node->cell.prob + i_node->cell.u[q]*r->cell.prob);
        }

        r->cell.dcu -= (G.dt/G.del[q])*(dcu_p-dcu_m);
    }
}

void BST::update_ctu(const Grid& G, TreeNode* r){
    if (r==nullptr){
        return;
    }
    update_ctu(G,r->left);
    update_ctu(G,r->right);

    for(int q =0; q < DIM; q++){
        TreeNode* i_node = r->cell.i_nodes[q];
        TreeNode* j_node; TreeNode* p_node;
        if(i_node!=dead){
            double F = G.dt*(r->cell.prob-i_node->cell.prob)/(2*G.del[q]);
                for(int b = 0; b < DIM; b++){
                    if (b!=q){
                        j_node = r->cell.i_nodes[b];
                        p_node = i_node->cell.i_nodes[b];

                        r->cell.ctu[b]      -= i_node->cell.w[q] * r->cell.w[b] * F;
                        j_node->cell.ctu[b] -= i_node->cell.w[q] * j_node->cell.u[b] * F;
                        i_node->cell.ctu[b] -= i_node->cell.u[q] * i_node->cell.w[b] * F;
                        p_node->cell.ctu[b] -= i_node->cell.u[q] * p_node->cell.u[b] * F;
                    }
                }

            // High-Resolution Correction Terms
            double th;
            if (i_node->cell.v[q]>0){
                TreeNode* i_i_node = i_node->cell.i_nodes[q];
                th = (i_node->cell.prob-i_i_node->cell.prob)/(r->cell.prob-i_node->cell.prob);
            }else{
                th = (r->cell.k_nodes[q]->cell.prob-r->cell.prob)/(r->cell.prob-i_node->cell.prob);
            }

            i_node->cell.ctu[q] += std::abs(i_node->cell.v[q])*(G.del[q]/G.dt - std::abs(i_node->cell.v[q]))*F*mc(th);
        }
    }
}

void BST::update_prob(const Grid& G, TreeNode* r){
    if (r == nullptr){
        return;
    }

    update_prob(G, r->left);
    update_prob(G, r->right);

    r->cell.prob += r->cell.dcu;
    for(int q =0; q < DIM; q++){
        r->cell.prob -= (G.dt/G.del[q])*(r->cell.ctu[q]-r->cell.i_nodes[q]->cell.ctu[q]);
    }
}

void BST::record_data(const std::string& file_name, const Grid& G, double t){
    std::ofstream myfile; myfile.open(file_name);
    myfile << t << std::endl;
    write_file(myfile, G, root);
    myfile.close();
}

void BST::max_key(const Grid& G, TreeNode* r, uint64_t& key){
    if (r == nullptr){
        return;
    }

    max_key(G, r->left, key);
    max_key(G, r->right, key);

    key = std::max(key, r->key);
}

void BST::measurement_update(Measurement m, const Grid& G, TreeNode* r){

    if (r == nullptr){
        return;
    }

    measurement_update(m, G, r->left);
    measurement_update(m, G, r->right);

    double x = 0;
    for(int q =0; q < DIM; q++){
        if (m.std[q] != 0){
            x += pow((r->cell.state[q]*G.del[q]-m.mean[q]),2)/pow(m.std[q],2);
        }
    }

    double prob = exp(-x/2); r->cell.prob *= prob;
}

void BST::check_cfl_condition(const Grid& G, TreeNode* r){
    if (r == nullptr){
        return;
    }

    check_cfl_condition(G, r->left);
    check_cfl_condition(G, r->right);

    cfl_min_dt = std::min(cfl_min_dt,r->cell.cfl_dt);
}

#endif //GBEES_BST_H
