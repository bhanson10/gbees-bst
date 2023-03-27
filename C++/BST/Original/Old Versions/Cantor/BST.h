#ifndef BST_h
#define BST_h

#include <iostream>
#include <wchar.h>
#include <array>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>
const int DIM = 3; 
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
class TreeNode; 

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
        std::array<TreeNode*, DIM> i_nodes; 
        std::array<TreeNode*, DIM> k_nodes; 
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

class TreeNode{ // TreeNode - Node in a Binary search tree
    public:
        int key;
        Cell cell; 
        TreeNode* left; 
        TreeNode* right; 
        TreeNode(int k, Cell c);
}; 

class BST{ // Binary Search Tree - data structure storing information
    public: 
        TreeNode* root; 

        BST();
        TreeNode* insertRecursive(TreeNode* r, TreeNode* new_node);
        TreeNode* recursiveSearch(TreeNode* r, int k);
        TreeNode* minValueNode(TreeNode* node);
        TreeNode* deleteNode(TreeNode* r, int k);
        int getHeight(TreeNode* r);
        int getDifference(TreeNode* r);
        TreeNode* rr_rotate(TreeNode* parent);
        TreeNode* ll_rotate(TreeNode* parent);
        TreeNode* lr_rotate(TreeNode* parent);
        TreeNode* rl_rotate(TreeNode* parent);
        TreeNode* balance(TreeNode* r);
        bool isBalanced(TreeNode* r);
        void writeFile(std::ofstream& myfile, Grid G, TreeNode* r);
        int full_size(TreeNode* r);
        int active_size(TreeNode* r, Grid G);
};

/*==============================================================================
FUNCTION DEFINITIONS
==============================================================================*/
TreeNode::TreeNode(int k, Cell c){ // Initialization of Node
    key = k; 
    cell = c;
    left = NULL; right = NULL; 
};

BST::BST(){ // Initialization of BST
    root = NULL; 
};

TreeNode* BST::insertRecursive(TreeNode* r, TreeNode* new_node){ // Inserting a node into BST
    if (r==NULL){
        r=new_node;
        return r; 
    }
    if (new_node->key < r->key){
        r->left = insertRecursive(r->left, new_node);
    }else if (new_node->key > r->key){
        r->right = insertRecursive(r->right, new_node);
    }else{
        r->cell = new_node->cell; 
        return r; 
    }
    return r; 
};

TreeNode* BST::recursiveSearch(TreeNode* r, int k){ // Searching to see if a key exists in BST
    if ((r==NULL)||(r->key == k)){
        return r;
    }else if (k < r->key){
        return recursiveSearch(r->left, k);
    }else{
        return recursiveSearch(r->right,k);
    }
};

TreeNode* BST::minValueNode(TreeNode* node){ // Getting the minimum value of a node for deletion
    TreeNode* current = node; 
    while(current->left != NULL){
        current = current->left; 
    }
    return current;
};

TreeNode* BST::deleteNode(TreeNode* r, int k){ // Deleting a node from the BST
    if(r==NULL){
        return NULL;
    }else if(k < r->key){
        r->left = deleteNode(r->left, k);
    }else if(k > r->key){
        r->right = deleteNode(r->right,k);
    }else{
        if(r->left == NULL){
            TreeNode* temp = r->right;
            delete r; 
            return temp; 
        } else if (r->right == NULL){
            TreeNode* temp = r->left; 
            delete r; 
            return temp; 
        }else{
            TreeNode* temp = minValueNode(r->right);
            r->key = temp->key; r->cell = temp->cell; 
            r->right = deleteNode(r->right, temp->key);
        }
    }
    return r;
};

int BST::getHeight(TreeNode* r){
    if (r == NULL){
        return 0;
    }

    return 1 + std::max(getHeight(r->left), getHeight(r->right)); 
};

int BST::getDifference(TreeNode* r){
    int l_height = getHeight(r->left);
    int r_height = getHeight(r->right);
    int b_factor = l_height - r_height;
    return b_factor; 
};

TreeNode* BST::rr_rotate(TreeNode* parent){
    TreeNode* t; 
    t = parent->right;
    parent->right = t->left;
    t->left = parent; 
    return t; 
};

TreeNode* BST::ll_rotate(TreeNode* parent){
    TreeNode* t; 
    t = parent->left;
    parent->left = t->right;
    t->right = parent; 
    return t; 
};

TreeNode* BST::lr_rotate(TreeNode* parent){
    TreeNode* t; 
    t = parent->left;
    parent->left = rr_rotate(t);
    return ll_rotate(parent); 
};

TreeNode* BST::rl_rotate(TreeNode* parent){
    TreeNode* t; 
    t = parent->right;
    parent->right = ll_rotate(t);
    return rr_rotate(parent); 
};

TreeNode* BST::balance(TreeNode* r){
    int bal_factor = getDifference(r); 
    while(abs(bal_factor) > 1){
        if (bal_factor > 1){
            if (getDifference(r->left) > 0){
                r = ll_rotate(r);
            }else{
                r = lr_rotate(r);
            }
        }else if (bal_factor < -1){
            if (getDifference(r->right) > 0){
                r = rl_rotate(r);
            }else{
                r = rr_rotate(r);
            }
        }
        bal_factor = getDifference(r);
    }
    return r; 
};

bool BST::isBalanced(TreeNode* r){
    int bal_factor = getDifference(r);
    if(abs(bal_factor) > 1){
        return 0;
    }else{
        return 1; 
    }
};

int BST::full_size(TreeNode* r){
    if(r == NULL){
        return 0;
    }else{
        return 1 + full_size(r->left) + full_size(r->right);
    }
};

int BST::active_size(TreeNode* r, Grid G){
    if(r == NULL){
        return 0;
    }else{
        if(r->cell.prob >= G.thresh){
            return 1 + active_size(r->left,G) + active_size(r->right,G);
        }else{
            return 0 + active_size(r->left,G) + active_size(r->right,G);
        }
    }
};

void BST::writeFile(std::ofstream& myfile, Grid G, TreeNode* r){
    if (r == NULL){
        return; 
    }

    writeFile(myfile, G, r->left);
    writeFile(myfile, G, r->right);

    myfile << r->cell.prob << " " << r->cell.state[0] << " " << r->cell.state[1] << " " << r->cell.state[2] << " " << r->cell.K << std::endl;
};

/*
void BST::writeFile(std::ofstream& myfile, Grid G, TreeNode* r){
    if (r == NULL){
        return; 
    }

    writeFile(myfile, G, r->left);
    writeFile(myfile, G, r->right);

    if (r->cell.prob >= G.thresh){
        myfile << r->cell.prob << " " << r->cell.state[0] << " " << r->cell.state[1] << " " << r->cell.state[2] << std::endl;
    }
};
*/
#endif // BST_h