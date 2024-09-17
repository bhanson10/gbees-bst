// gbees.h, https://github.com/bhanson10/gbees
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef GBEES_H
#define GBEES_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#define TOL 1E-8

/*==============================================================================
                            STRUCTURE DEFINITIONS
==============================================================================*/
typedef struct Meas {
    int dim; 
    double *mean;
    double **cov;
    double T;
} Meas;

Meas Meas_create(int dim, const char* M_DIR, const char* M_FILE);

void Meas_free(Meas *M);

typedef struct Grid {
    int dim; 
    double thresh;
    double dt;
    double *center;
    double *dx;
    double hi_bound;
    double lo_bound;
} Grid;

Grid Grid_create(int dim, double thresh, double* center, double* dx);

void Grid_free(Grid* G);

typedef struct Traj {
    double *coef;
} Traj;

Traj Traj_create(int n, double* coef);

void Traj_free(Traj* T);

typedef struct TreeNode TreeNode;

struct TreeNode { // REF- Remove duplicated typedef
    uint64_t key;
    double prob;
    double *v;
    double *ctu;
    int *state;
    TreeNode **i_nodes;
    TreeNode **k_nodes;
    double dcu;
    double cfl_dt;
    int new_f;
    int ik_f;
    double bound_val; 
    TreeNode* left;
    TreeNode* right;
};

TreeNode* TreeNode_create(int dim, uint64_t key, double prob, int* state, double J);

void Tree_free(TreeNode* r);

void TreeNode_free(TreeNode* r);



/*==============================================================================
                        NON-MEMBER FUNCTION DEFINITIONS
==============================================================================*/

uint64_t rosenberg_pair(const int* state, int d, int m);

uint64_t state_conversion(int dim, const int* state);

double mc(double th);

void inv_mat(double** mat, double* inverse, int size);

void mul_mat_vec(double* matrix, double* vector, double* result, int size);

double dot_product(double* vec1, double* vec2, int size);

double gauss_probability(int dim, double* x, Meas M);

/*==============================================================================
                    BINARY SEARCH TREE FUNCTION DEFINITIONS
==============================================================================*/
TreeNode* insert_recursive(TreeNode* r, Grid* G, Traj T, double prob, uint64_t key, int* state, bool BOUNDS, double (*BOUND_f)(double*, double*));

TreeNode* search_recursive(TreeNode* r, uint64_t key);

TreeNode* min_value_node(TreeNode* node);

TreeNode* delete_node(TreeNode* r, uint64_t k, Grid G);

int get_height(TreeNode* r);

int get_difference(TreeNode* r);

TreeNode* rr_rotate(TreeNode* parent);

TreeNode* ll_rotate(TreeNode* parent);

TreeNode* rl_rotate(TreeNode* parent);

TreeNode* lr_rotate(TreeNode* parent);

TreeNode* balance(TreeNode* r);


/*==============================================================================
                        GBEES FUNCTION DEFINITIONS
==============================================================================*/
void initialize_adv(void (*f)(double*, double*, double*, double*), TreeNode* r, Grid* G, Traj T);

void initialize_ik_nodes(TreeNode* P, TreeNode* r, Grid* G);

void recursive_loop(TreeNode** P, Grid* G, Meas M, Traj T, int level, int* current_state, double* current_state_vec, bool BOUNDS, double (*BOUND_f)(double*, double*));

void initialize_grid(void (*f)(double*, double*, double*, double*), TreeNode** P, Grid* G, Meas M, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*));

void set_bounds(TreeNode* r, Grid* G);

void get_sum(TreeNode* r, double* prob_sum);

void divide_sum(TreeNode* r, double prob_sum, Grid* G, uint64_t* max_key, int* a_count, int* tot_count);

void normalize_tree(TreeNode* P, Grid* G, uint64_t* max_key, int* a_count, int* tot_count);

char* concat_m(const char* str1, const char* str2, int num1);

char* concat_p(const char* str1, const char* str2, int num1, const char* str3, int num2);

void write_file(FILE* myfile, TreeNode* r, Grid G);

void record_data(TreeNode* r, const char* FILE_NAME, Grid G, const double t);

void create_neighbors(TreeNode** P, TreeNode* r, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*));

void grow_tree(void (*f)(double*, double*, double*, double*), TreeNode** P, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*));

void check_cfl_condition(TreeNode* r, Grid* G);

void update_dcu(TreeNode* r, Grid G);

void update_ctu(TreeNode* r, Grid G);

void godunov_method(TreeNode** P, Grid G);

void update_prob(TreeNode* r, Grid G);

void mark_cells(TreeNode* r, Grid G, double* del_probs, uint64_t* del_keys, int* idx);

// REF- different definition of qsort_r in MAC and Linux
#ifdef __linux__ 
  int compare_indices(const void *a, const void *b, void *del_probs);
#else  
  int compare_indices(void *del_probs, const void *a, const void *b);
#endif  

void sort_by_double(double *del_probs, uint64_t *del_keys, size_t n);

void delete_cells(TreeNode** P, Grid G, double* del_probs, uint64_t* del_keys, int idx);

void prune_tree(TreeNode** P, Grid G, int tot_count);

void meas_up_recursive(void (*h)(double*, double*, double*, double*), TreeNode* r, Grid G, Meas M, Traj T);

void run_gbees(void (*f)(double*, double*, double*, double*), void (*h)(double*, double*, double*, double*), double (*BOUND_f)(double*, double*), Grid G, Meas M, Traj T, char* P_DIR, char* M_DIR, int NUM_DIST, int NUM_MEAS, int DEL_STEP, int OUTPUT_FREQ, int DIM_h, bool OUTPUT, bool RECORD, bool MEASURE, bool BOUNDS);

#endif
