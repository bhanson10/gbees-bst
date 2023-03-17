#include "header.h"
#include <ctime>
#include <cstdlib>

int main()
{   
    BST P; double prob; double prob_sum = 0; int key; 
    srand((unsigned)time(NULL)); int n = 5; 

    for (int i = 0; i < n; i++){
        //key = 1 + (rand()%100);
        key = i; 

        prob = (float)rand()/RAND_MAX; prob_sum += prob; 
        Cell c = {.prob=prob,.active=1};
        TreeNode* new_node = new TreeNode(key, c);
        P.root = P.insertRecursive(P.root, new_node);
    }
    P.normalizeProb(P.root, prob_sum);

    TreeNode* test = P.recursiveSearch(P.root, 2);
    cout << test->key << endl; 
    cout << test->cell.prob << endl;
    P.print2D(P.root, 1);

    P.root = P.balance(P.root);
    cout << test->key << endl; 
    cout << test->cell.prob << endl;
    P.print2D(P.root, 1);
    
    /*
    cout << "--------------------------------------------" << endl; 
    P.print2D(P.root, 1);
    cout << endl; 
    cout << "HEIGHT: " << P.getHeight(P.root) << endl;
    cout << "--------------------------------------------" << endl; 

    TreeNode* test = P.recursiveSearch(P.root, 2);
    cout << test->key << endl; 
    cout << test->cell.prob << endl;

    P.root = P.delete_single(P.root, test);
    cout << "--------------------------------------------" << endl; 
    P.print2D(P.root, 1);
    cout << endl; 
    cout << "HEIGHT: " << P.getHeight(P.root) << endl;
    cout << "--------------------------------------------" << endl; 

    cout << test->key << endl;
    cout << test->cell.prob << endl; 
    */
    
    /*
    if(!P.isBalanced(P.root)){
        cout << "NOT BALANCED... BALANCING" << endl; 
        P.root = P.balance(P.root);
        cout << "--------------------------------------------" << endl; 
        P.print2D(P.root, 1);


        cout << endl; 
        cout << "HEIGHT: " << P.getHeight(P.root) << endl;
        cout << "--------------------------------------------" << endl; 
        cout << "BALANCED" << endl; 
    }else{
        cout << "BALANCED" << endl; 
    }
    
    P.root = P.addNodes(P.root); 
    cout << "ADDING NODES WITH PROBABILITIES > 0.1..." << endl; 
    cout << "--------------------------------------------" << endl; 
    P.print2D(P.root, 1);
    cout << endl; 
    cout << "HEIGHT: " << P.getHeight(P.root) << endl;
    cout << "--------------------------------------------" << endl; 
    
    if(!P.isBalanced(P.root)){
        cout << "NOT BALANCED... BALANCING" << endl; 
        P.root = P.balance(P.root);
        cout << "--------------------------------------------" << endl; 
        P.print2D(P.root, 1);
        cout << endl; 
        cout << "HEIGHT: " << P.getHeight(P.root) << endl;
        cout << "--------------------------------------------" << endl; 
        cout << "BALANCED" << endl; 
    }else{
        cout << "BALANCED" << endl; 
    }
    P.root = P.deleteNodes(P.root);
    cout << "DELETING NODES WITH PROBABILITIES < 0.1 (not including nodes just added)..." << endl; 
    cout << "--------------------------------------------" << endl; 
    P.print2D(P.root, 1);
    cout << endl; 
    cout << "HEIGHT: " << P.getHeight(P.root) << endl;
    cout << "--------------------------------------------" << endl; 

    if(!P.isBalanced(P.root)){
        cout << "NOT BALANCED... BALANCING" << endl; 
        P.root = P.balance(P.root);
        cout << "--------------------------------------------" << endl; 
        P.print2D(P.root, 1);
        cout << endl; 
        cout << "HEIGHT: " << P.getHeight(P.root) << endl;
        cout << "--------------------------------------------" << endl; 
        cout << "BALANCED" << endl; 
    }else{
        cout << "BALANCED" << endl; 
    }
    */
    return 0;
}