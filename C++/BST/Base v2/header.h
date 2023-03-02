#ifndef header_h
#define header_h

#include <iostream>
#include <array>
#include <cmath> 
using namespace std; 
#define SPACE 10

const int DIM = 3; 
const int DEAD = -1; 

struct Cell{  // Cell Class
    double prob;
    int active; 
};

class TreeNode{
    public:
        int key;
        Cell cell; 
        TreeNode* left; 
        TreeNode* right; 

        TreeNode(int k, Cell c){
            key = k; 
            cell = {.prob=c.prob,.active=c.active};
            left = NULL; right = NULL; 
        }
}; 

class BST{
    public: 
        TreeNode* root; 
        TreeNode* dead = new TreeNode(-1, {.prob = 0, .active = 0}); 

        BST(){ // Initialization of BST
            root = NULL; 
        }

        TreeNode* insertRecursive(TreeNode* r, TreeNode* new_node){ // Recursively inserting a node into BST
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
        }

        TreeNode* recursiveSearch(TreeNode* r, int k){ // Searching to see if a key exists in BST
            if ((r==NULL)||(r->key == k)){
                return r;
            }else if (k < r->key){
                return recursiveSearch(r->left, k);
            }else{
                return recursiveSearch(r->right,k);
            }
        }

        TreeNode* minValueNode(TreeNode* node){ // Getting the minimum value of a node for deletion
            TreeNode* current = node; 
            while(current->left != NULL){
                current = current->left; 
            }
            return current;
        }

        TreeNode* delete_single(TreeNode* r, int k){ // Deleting a node from the BST
            if(r==NULL){
                return NULL;
            }else if(k < r->key){
                r->left = delete_single(r->left, k);
            }else if(k > r->key){
                r->right = delete_single(r->right,k);
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
                    r->right = delete_single(r->right, temp->key);
                }
            }
            return r;
        }

        void normalizeProb(TreeNode* r, double prob_sum){
            if (r==NULL){
                return;
            }
            normalizeProb(r->left, prob_sum);
            r->cell.prob /= prob_sum; //Perform the Action
            normalizeProb(r->right, prob_sum); 
        };

        void print2D(TreeNode* r, int space){ // Printing out all the nodes in a BST in a 2D fashion
            if (r == NULL){
                return;
            }
            space+= SPACE;
            print2D(r->right, space);
            cout << endl; 
            for(int i = SPACE; i < space; i++){
                cout << " ";
            }
            cout << r->key << ", " << r->cell.prob << "\n";
            print2D(r->left, space);
        };

        void printInorder(TreeNode* r){ //  (Left, current node, Right)
            if (r == NULL){
                return;
            }
            printInorder(r->left);
            std::cout << r->key << " " << r->cell.prob << std::endl;
            printInorder(r -> right);
        };

        TreeNode* addNodes(TreeNode* r){
            if(r==NULL){
                return r; 
            }
            addNodes(r->left);
            addNodes(r->right);

            if(r->cell.prob > 0.1){
                int key = r->key + 100; 
                Cell c = {.prob=0,.active=0};
                TreeNode* new_node = new TreeNode(key, c);
                root = insertRecursive(root, new_node);
            }
            return root; 
        };

        TreeNode* deleteNodes(TreeNode* r){
            if(r==NULL){
                return r; 
            }
            deleteNodes(r->left);
            deleteNodes(r->right);

            if((r->cell.prob < 0.1)&&(r->cell.active == 1)){
                root = delete_single(root, r->key);
            }
            return root; 
        };

        int getHeight(TreeNode* r){
            if (r == NULL){
                return 0;
            }

            return 1 + max(getHeight(r->left), getHeight(r->right)); 
        };

        int getDifference(TreeNode* r){
            int l_height = getHeight(r->left);
            int r_height = getHeight(r->right);
            int b_factor = l_height - r_height;
            return b_factor; 
        };

        TreeNode* rr_rotate(TreeNode* parent){
            TreeNode* t; 
            t = parent->right;
            parent->right = t->left;
            t->left = parent; 
            return t; 
        };

        TreeNode* ll_rotate(TreeNode* parent){
            TreeNode* t; 
            t = parent->left;
            parent->left = t->right;
            t->right = parent; 
            return t; 
        };

        TreeNode* lr_rotate(TreeNode* parent){
            TreeNode* t; 
            t = parent->left;
            parent->left = rr_rotate(t);
            return ll_rotate(parent); 
        };

        TreeNode* rl_rotate(TreeNode* parent){
            TreeNode* t; 
            t = parent->right;
            parent->right = ll_rotate(t);
            return rr_rotate(parent); 
        };

        TreeNode* balance(TreeNode* r){
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

        bool isBalanced(TreeNode* r){
            int bal_factor = getDifference(r);
            if(abs(bal_factor) > 1){
                return 0;
            }else{
                return 1; 
            }
        };
};

#endif // header_h