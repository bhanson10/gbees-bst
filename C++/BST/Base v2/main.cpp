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

    P.print2D(P.root, 1); 

    P.delete_single(P.root,3);
    P.print2D(P.root, 1); 
    P.delete_single(P.root,6);
    P.print2D(P.root, 1); 
    return 0;
}