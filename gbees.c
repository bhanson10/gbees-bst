// gbees.c, https://github.com/bhanson10/gbees
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <time.h>

#include "gbees.h"


Meas Meas_create(int dim, const char* M_DIR, const char* M_FILE) {
    char M_PATH[256];
    snprintf(M_PATH, sizeof(M_PATH), "%s/%s", M_DIR, M_FILE);

    FILE *m_file = fopen(M_PATH, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s", M_PATH);
        exit(EXIT_FAILURE);
    }

    Meas M;
    M.dim = dim; 
    M.mean = malloc(dim * sizeof(double));
    M.cov = malloc(dim * sizeof(double *));
    if (M.mean == NULL || M.cov == NULL) {
        perror("Error: memory allocation failure during measurement creation");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < dim; i++) {
        M.cov[i] = malloc(dim * sizeof(double));
        if (M.cov[i] == NULL) {
            perror("Error: memory allocation failure during measurement creation");
            exit(EXIT_FAILURE);
        }
    }

    char line[256];
    char *token;
    int count = 0;

    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // mean vector
    token = strtok(line, " ");
    while (token != NULL && count < dim) { // read mean vector
        M.mean[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0;

    // Read covariance matrix
    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    for (int i = 0; i < dim; i++) { // read covariance matrix
        fgets(line, sizeof(line), m_file);
        token = strtok(line, " ");
        while (token != NULL && count < dim) {
            M.cov[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0;
    }

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // read T value
    M.T = strtod(line, NULL);

    fclose(m_file);
    return M;
}

void Meas_free(Meas *M) {
    if (M->mean != NULL) {
        free(M->mean);
    }
    if (M->cov != NULL) {
        for (int i = 0; i < M->dim; i++) {
            if (M->cov[i] != NULL) {
                free(M->cov[i]);
            }
        }
        free(M->cov);
    }
}

Grid Grid_create(int dim, double thresh, double* center, double* dx){
    Grid G; 
    G.dim = dim; 
    G.thresh = thresh; 
    G.dt = DBL_MAX; 
    G.center = malloc(dim * sizeof(double));
    G.dx = malloc(dim * sizeof(double *));
    if (G.center == NULL || G.dx == NULL) {
        perror("Error: memory allocation failure during grid creation");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < dim; i++){
        G.center[i] = center[i]; 
        G.dx[i] = dx[i]; 
    }
    G.hi_bound = DBL_MAX; 
    G.lo_bound = -DBL_MAX; 
    return G;
}

void Grid_free(Grid* G) {
    // Free the allocated memory
    free(G->center);
    free(G->dx);

    // Set pointers to NULL to avoid dangling pointers
    G->center = NULL;
    G->dx = NULL;
}


Traj Traj_create(int n, double* coef){
    Traj T; 
    T.coef = malloc(n * sizeof(double));
    if (T.coef == NULL) {
        perror("Error: memory allocation failure during trajectory creation");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < n; i++){
        T.coef[i] = coef[i]; 
    }
    return T;
}

void Traj_free(Traj* T) {
    if (T->coef != NULL) {
        free(T->coef);
        T->coef = NULL;
    }
}

TreeNode* TreeNode_create(int dim, uint64_t key, double prob, int* state, double J) {
    TreeNode* node = (TreeNode*)malloc(sizeof(TreeNode));
    if (node == NULL) {
        perror("Error: memory allocation failure during TreeNode creation");
        exit(EXIT_FAILURE);
    }
    node->key = key;
    node->prob = prob; 
    node->v = (double*)malloc(dim * sizeof(double));
    node->ctu = (double*)malloc(dim * sizeof(double));
    node->state = (int*)malloc(dim * sizeof(int));
    node->i_nodes = (TreeNode **)malloc(dim * sizeof(TreeNode *));
    node->k_nodes = (TreeNode **)malloc(dim * sizeof(TreeNode *));

    // Check for memory allocation failure
    if (node->v == NULL || node->ctu == NULL || node->state == NULL || node->i_nodes == NULL || node->k_nodes == NULL) {
        if (node->v) free(node->v);
        if (node->ctu) free(node->ctu);
        if (node->state) free(node->state);
        if (node->i_nodes) free(node->i_nodes);
        if (node->k_nodes) free(node->k_nodes);
        perror("Error: memory allocation failure during TreeNode creation");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < dim; i++) {
        node->v[i] = 0.0;
        node->ctu[i] = 0.0;
        node->state[i] = state[i];
        node->i_nodes[i] = NULL;
        node->k_nodes[i] = NULL;
    }
    node->new_f = 0; 
    node->bound_val = J; 
    node->left = NULL; 
    node->right = NULL; 
    return node;
}

void Tree_free(TreeNode* r) {
    if (r == NULL) {
        return;
    }
    
    Tree_free(r->left);
    Tree_free(r->right);
    free(r->v);
    free(r->ctu);
    free(r->state);
    free(r->i_nodes);
    free(r->k_nodes);
    free(r);
}

void TreeNode_free(TreeNode* r) {
    free(r->v);
    free(r->ctu);
    free(r->state);
    free(r->i_nodes);
    free(r->k_nodes);
    free(r);
}

/*==============================================================================
                        NON-MEMBER FUNCTION DEFINITIONS
==============================================================================*/
uint64_t rosenberg_pair(const int* state, int d, int m) { // Recursive Rosenberg pairing
    if (d == 1) {
        return state[0];
    }

    int* new_state = (int*)malloc((d-1) * sizeof(int));
    if (new_state == NULL) {
        perror("Error: memory allocation failure during Rosenberg pairing");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < d-1; i++) {
        new_state[i] = state[i];
    }

    int new_m = new_state[0];
    for (int i = 1; i < d-1; i++) {
        if (new_state[i] > new_m) {
            new_m = new_state[i];
        }
    }

    uint64_t result = rosenberg_pair(new_state, d-1, new_m) + (uint64_t)pow(m, d) + (m - state[d-1]) * ((uint64_t)pow(m+1, d-1) - (uint64_t)pow(m, d-1));
    free(new_state);
    return result;
}

uint64_t state_conversion(int dim, const int* state) { // Convert from n-dimensional state to unique key via Rosenberg pairing
    int* shift_state = (int*)malloc(dim * sizeof(int));
    if (shift_state == NULL) {
        perror("Error: memory allocation failure during state conversion");
        exit(EXIT_FAILURE);
    }

    int m;
    uint64_t key;

    // Perform negative shift
    for (int i = 0; i < dim; i++) {
        if (state[i] < 0) {
            shift_state[i] = -2 * state[i] - 1;
        } else {
            shift_state[i] = 2 * state[i];
        }
    }

    // Determine the maximum value in shift_state
    m = shift_state[0];
    for (int i = 1; i < dim; i++) {
        if (shift_state[i] > m) {
            m = shift_state[i];
        }
    }

    // Compute the key using the Rosenberg pairing function
    key = rosenberg_pair(shift_state, dim, m);

    // Free the allocated memory
    free(shift_state);
    
    return key;
}

double mc(double th){ // MC flux limiter
    double min1 = fmin((1 + th)/2.0, 2.0);
    return fmax(0.0, fmin(min1, 2*th)); 
}

void inv_mat(double** mat, double* inverse, int size) {
    int i, j, k;
    double ratio;
    double* augmented = (double*)malloc(size * size * 2 * sizeof(double));

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            augmented[i * 2 * size + j] = mat[i][j];
            augmented[i * 2 * size + (j + size)] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (i = 0; i < size; i++) {
        if (augmented[i * 2 * size + i] == 0) {
            perror("Error: matrix inversion error, zero pivot element");
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < size; j++) {
            if (i != j) {
                ratio = augmented[j * 2 * size + i] / augmented[i * 2 * size + i];
                for (k = 0; k < 2 * size; k++) {
                    augmented[j * 2 * size + k] -= ratio * augmented[i * 2 * size + k];
                }
            }
        }
    }

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            inverse[i * size + j] = augmented[i * 2 * size + (j + size)] / augmented[i * 2 * size + i];
        }
    }

    free(augmented);
}

void mul_mat_vec(double* matrix, double* vector, double* result, int size) {
    int i, j;
    for (i = 0; i < size; i++) {
        result[i] = 0;
        for (j = 0; j < size; j++) {
            result[i] += matrix[i * size + j] * vector[j];
        }
    }
}

double dot_product(double* vec1, double* vec2, int size) {
    int i;
    double result = 0;
    for (i = 0; i < size; i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

double gauss_probability(int dim, double* x, Meas M){ // MC Calculate gaussian probability at state x given mean and covariance
    double* M_inv = (double*)malloc(dim * dim * sizeof(double));
    double* M_inv_x = (double*)malloc(dim * sizeof(double));
    double exp_result;
    
    double diff[dim];
    for(int i = 0; i < dim; i++){
        diff[i] = x[i] - M.mean[i]; 
    }
    inv_mat(M.cov, M_inv, dim);
    mul_mat_vec(M_inv, diff, M_inv_x, dim);
    double dot_prod = dot_product(diff, M_inv_x, dim);
    exp_result = exp(-0.5 * dot_prod);

    free(M_inv);
    free(M_inv_x);

    return exp_result;
}

/*==============================================================================
                    BINARY SEARCH TREE FUNCTION DEFINITIONS
==============================================================================*/
TreeNode* insert_recursive(TreeNode* r, Grid* G, Traj T, double prob, uint64_t key, int* state, bool BOUNDS, double (*BOUND_f)(double*, double*)) {
    if (r == NULL) {
        if(BOUNDS){
            double x[G->dim];
            for(int k = 0; k < G->dim; k++){
                x[k] = G->dx[k]*state[k]+G->center[k];
            }

            double J = BOUND_f(x, T.coef); 

            if(J >= G->lo_bound && J <= G->hi_bound){
                return TreeNode_create(G->dim, key, prob, state, J);
            }else{
                return NULL; 
            }
        }
        else{
            return TreeNode_create(G->dim, key, prob, state, 0.0);
        }
    }

    if (key < r->key) {
        r->left = insert_recursive(r->left, G, T, prob, key, state, BOUNDS, BOUND_f);
    } else if (key > r->key) {
        r->right = insert_recursive(r->right, G, T, prob, key, state, BOUNDS, BOUND_f);
    }

    return r;
}

TreeNode* search_recursive(TreeNode* r, uint64_t key) {
    if (r == NULL || r->key == key) {
        return r;
    }

    if (key < r->key) {
        return search_recursive(r->left, key);
    } else {
        return search_recursive(r->right, key);
    }
}


TreeNode* min_value_node(TreeNode* node){ 
    TreeNode* current = node;
    while(current->left != NULL){
        current = current->left;
    }
    return current;
}

TreeNode* delete_node(TreeNode* r, uint64_t k, Grid G){ 
    if(r==NULL){
        return NULL;
    }else if(k < r->key){
        r->left = delete_node(r->left, k, G);
    }else if(k > r->key){
        r->right = delete_node(r->right,k, G);
    }else{
        if(r->left == NULL){
            TreeNode* temp = r->right;
            TreeNode_free(r);
            return temp;
        } else if (r->right == NULL){
            TreeNode* temp = r->left;
            TreeNode_free(r); 
            return temp;
        }else{
            TreeNode* temp = min_value_node(r->right);
            r->key = temp->key; 
            r->prob = temp->prob; 
            for(int i = 0; i < G.dim; i++){
                r->v[i] = temp->v[i]; 
                r->ctu[i] = temp->ctu[i]; 
                r->state[i] = temp->state[i]; 
                r->i_nodes[i] = temp->i_nodes[i]; 
                r->k_nodes[i] = temp->k_nodes[i]; 
            }
            r->dcu = temp->dcu;
            r->cfl_dt = temp->cfl_dt;
            r->new_f = temp->new_f;
            r->ik_f = temp->ik_f;
            r->bound_val = temp->bound_val; 
            r->right = delete_node(r->right, temp->key, G);
        }
    }
    return r;
}

int get_height(TreeNode* r){
    if (r == NULL){
        return 0;
    }

    return 1 + fmax(get_height(r->left), get_height(r->right));
}

int get_difference(TreeNode* r){
    int l_height = get_height(r->left);
    int r_height = get_height(r->right);
    int b_factor = l_height - r_height;
    return b_factor;
}

TreeNode* rr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = t->left;
    t->left = parent;
    return t;
}

TreeNode* ll_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = t->right;
    t->right = parent;
    return t;
}

TreeNode* rl_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = ll_rotate(t);
    return rr_rotate(parent);
}

TreeNode* lr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = rr_rotate(t);
    return ll_rotate(parent);
}

TreeNode* balance(TreeNode* r){
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

/*==============================================================================
                        GBEES FUNCTION DEFINITIONS
==============================================================================*/
void initialize_adv(void (*f)(double*, double*, double*, double*), TreeNode* r, Grid* G, Traj T){
    if(r == NULL){ 
        return;
    }
    initialize_adv(f, r->left, G, T);
    initialize_adv(f, r->right, G, T);
    
    if(r->new_f==0){
        double x[G->dim];
        for(int i = 0; i < G->dim; i++){
            x[i] = G->dx[i]*r->state[i]+G->center[i];
        }
        double xk[G->dim];
        (*f)(xk, x, G->dx, T.coef); 

        double sum = 0;
        for(int i = 0; i < G->dim; i++){
            r->v[i] = xk[i];
            sum += fabs(r->v[i]) / G->dx[i];
        }
        r->new_f = 1;
        r->cfl_dt = 1.0/sum;
    }
}

void initialize_ik_nodes(TreeNode* P, TreeNode* r, Grid* G){
    if(r == NULL){
        return;
    }
    initialize_ik_nodes(P, r->left, G); 
    initialize_ik_nodes(P, r->right, G);

    if(r->ik_f == 0){
        for(int i = 0; i < G->dim; i++){
            // Initializing i, k nodes
            int i_state[G->dim]; memcpy(i_state, r->state, G->dim * sizeof(int)); i_state[i] = i_state[i] - 1; uint64_t i_key = state_conversion(G->dim, i_state); 
            int k_state[G->dim]; memcpy(k_state, r->state, G->dim * sizeof(int)); k_state[i] = k_state[i] + 1; uint64_t k_key = state_conversion(G->dim, k_state); 
            TreeNode* i_node = search_recursive(P, i_key); TreeNode* k_node = search_recursive(P, k_key); 
            r->i_nodes[i] = i_node; r->k_nodes[i] = k_node;

            if(i_node != NULL){
                i_node->k_nodes[i] = r; 
            }
            if(k_node != NULL){
                k_node->i_nodes[i] = r; 
            }
        }
        r->ik_f = 1; 
    } 
}


void recursive_loop(TreeNode** P, Grid* G, Meas M, Traj T, int level, int* current_state, double* current_state_vec, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    if (level == G->dim) {
        uint64_t key = state_conversion(G->dim, current_state);
        double prob = gauss_probability(M.dim, current_state_vec, M);
        *P = insert_recursive(*P, G, T, prob, key, current_state, BOUNDS, BOUND_f);
        return;
    }

    int start = (int) round(-3 * pow(M.cov[level][level], 0.5) / G->dx[level]);
    int end = (int) round(3 * pow(M.cov[level][level], 0.5) / G->dx[level]);
    for (int i = start; i <= end; i++) {
        current_state[level] = i;
        current_state_vec[level] = i * G->dx[level] + G->center[level];
        recursive_loop(P, G, M, T, level + 1, current_state, current_state_vec, BOUNDS, BOUND_f);
    }
}

void initialize_grid(void (*f)(double*, double*, double*, double*), TreeNode** P, Grid* G, Meas M, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    int current_state[G->dim]; double current_state_vec[G->dim];
    recursive_loop(P, G, M, T, 0, current_state, current_state_vec, BOUNDS, BOUND_f);
    *P = balance(*P); 
    initialize_adv(f, *P, G, T);
    initialize_ik_nodes(*P, *P, G);
}

void set_bounds(TreeNode* r, Grid* G){
    if(r == NULL){
        return;
    }
    set_bounds(r->left, G); 
    set_bounds(r->right, G);

    G->lo_bound = fmin(G->lo_bound, r->bound_val); 
    G->hi_bound = fmax(G->hi_bound, r->bound_val); 
}

void get_sum(TreeNode* r, double* prob_sum){
    if (r == NULL){
        return;
    }
    get_sum(r->left, prob_sum);
    get_sum(r->right, prob_sum);

    *prob_sum += r->prob;
}

void divide_sum(TreeNode* r, double prob_sum, Grid* G, uint64_t* max_key, int* a_count, int* tot_count){
    if (r == NULL){
        return;
    }
    divide_sum(r->left, prob_sum, G, max_key, a_count, tot_count);
    divide_sum(r->right, prob_sum, G, max_key, a_count, tot_count);

    r->prob /= prob_sum;

    *max_key = fmax(*max_key, r->key);
    if(r->prob >= G->thresh){
        *a_count += 1; 
    }
    *tot_count += 1; 
}

void normalize_tree(TreeNode* P, Grid* G, uint64_t* max_key, int* a_count, int* tot_count){
    *max_key = 0; *a_count = 0; *tot_count = 0; 
    double prob_sum = 0;
    get_sum(P, &prob_sum);
    divide_sum(P, prob_sum, G, max_key, a_count, tot_count);
}

char* concat_m(const char* str1, const char* str2, int num1) {
    int num1_len = snprintf(NULL, 0, "%d", num1);
    char* str3 = ".txt"; 
    size_t total_len = strlen(str1) + strlen(str2) + num1_len + strlen(str3) + 1; 
    char* result = (char*)malloc(total_len * sizeof(char));
    if (result == NULL) {
        perror("Error: memory allocation failure during concat_m");
        exit(EXIT_FAILURE);
    }
    strcpy(result, str1);
    strcat(result, str2);
    char num1_str[num1_len + 1];
    sprintf(num1_str, "%d", num1);
    strcat(result, num1_str);
    strcat(result, str3);
    
    return result;
}

char* concat_p(const char* str1, const char* str2, int num1, const char* str3, int num2) {
    int num1_len = snprintf(NULL, 0, "%d", num1);
    int num2_len = snprintf(NULL, 0, "%d", num2);
    char* str4 = ".txt"; 
    size_t total_len = strlen(str1) + strlen(str2) + num1_len + strlen(str3) + num2_len + strlen(str4) + 1; 
    char* result = (char*)malloc(total_len * sizeof(char));
    if (result == NULL) {
        perror("Error: memory allocation failure during concat_p");
        exit(EXIT_FAILURE);
    }
    strcpy(result, str1);
    strcat(result, str2);
    char num1_str[num1_len + 1];
    sprintf(num1_str, "%d", num1);
    strcat(result, num1_str);
    strcat(result, str3);
    char num2_str[num2_len + 1];
    sprintf(num2_str, "%d", num2);
    strcat(result, num2_str);
    strcat(result, str4);

    return result;
}

void write_file(FILE* myfile, TreeNode* r, Grid G){
    if (r == NULL){
        return;
    }

    write_file(myfile, r->left, G);
    write_file(myfile, r->right, G);

    if (r->prob >= G.thresh) {
        fprintf(myfile, "%.10e", r->prob);
        for (int i = 0; i < G.dim; i++) {
            fprintf(myfile, " %.10e", G.dx[i] * r->state[i] + G.center[i]);
        }
        fprintf(myfile, "\n");
    }
}

void record_data(TreeNode* r, const char* FILE_NAME, Grid G, const double t){
    FILE* file = fopen(FILE_NAME, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: could not open file %s", FILE_NAME);
        exit(EXIT_FAILURE);
    }

    fprintf(file, "%lf\n", t);
    write_file(file, r, G);
    fclose(file);
}

void create_neighbors(TreeNode** P, TreeNode* r, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    if (r == NULL){
        return;
    }
    create_neighbors(P, r->left, G, T, BOUNDS, BOUND_f);
    create_neighbors(P, r->right, G, T, BOUNDS, BOUND_f);
    
    if (r->prob >= G.thresh){
        double current_v[G.dim];
        memcpy(current_v, r->v, G.dim * sizeof(double));
        int current_state[G.dim]; 
        memcpy(current_state, r->state, G.dim * sizeof(int));
        int new_state[G.dim]; uint64_t new_key;
        for (int i = 0; i < G.dim; i++){
            memcpy(new_state, current_state, G.dim * sizeof(int));
            // Checking Forward Faces
            if(current_v[i] > 0){
                if(r->k_nodes[i] == NULL){
                    new_state[i] += 1;
                    new_key = state_conversion(G.dim, new_state);
                    *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] += 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                new_state[j] += 1;
                                new_key = state_conversion(G.dim, new_state);
                                *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                            }else if (current_v[j] < 0){
                                new_state[j] -= 1;
                                new_key = state_conversion(G.dim, new_state);
                                *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] += 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                if(r->k_nodes[i]->k_nodes[j] == NULL){
                                    new_state[j] += 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                                }
                            }else if (current_v[j] < 0){
                                if(r->k_nodes[i]->i_nodes[j] == NULL){
                                    new_state[j] -= 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                                }
                            }
                        }
                    }
                }
                // Checking Backward Faces
            }else if(current_v[i] < 0){
                if(r->i_nodes[i] == NULL){
                    new_state[i] -= 1;
                    new_key = state_conversion(G.dim, new_state);
                    *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] -= 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                new_state[j] += 1;
                                new_key = state_conversion(G.dim, new_state);
                                *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f); 

                            }else if (current_v[j] < 0){
                                new_state[j] -= 1;
                                new_key = state_conversion(G.dim, new_state);
                                *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] -= 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                if(r->i_nodes[i]->k_nodes[j] == NULL){
                                    new_state[j] += 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                                }
                            }else if (current_v[j] < 0){
                                if(r->i_nodes[i]->i_nodes[j] == NULL){
                                    new_state[j] -= 1;
                                    new_key = state_conversion(G.dim, new_state); 
                                    *P = insert_recursive(*P, &G, T, 0, new_key, new_state, BOUNDS, BOUND_f);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void grow_tree(void (*f)(double*, double*, double*, double*), TreeNode** P, Grid G, Traj T, bool BOUNDS, double (*BOUND_f)(double*, double*)){
    create_neighbors(P, *P, G, T, BOUNDS, BOUND_f);
    *P = balance(*P); 
    initialize_adv(f, *P, &G, T);
    initialize_ik_nodes(*P, *P, &G);
}

void check_cfl_condition(TreeNode* r, Grid* G){
    if (r == NULL){
        return;
    }

    check_cfl_condition(r->left, G);
    check_cfl_condition(r->right, G);

    G->dt = fmin(G->dt,r->cfl_dt);
}

void update_dcu(TreeNode* r, Grid G){
    if (r == NULL){
        return;
    }
    update_dcu(r->left, G);
    update_dcu(r->right, G);

    r->dcu = 0; 
    TreeNode* i_node; TreeNode* k_node; 
    for(int i = 0; i < G.dim; i++){
        r->ctu[i] = 0.0; 
        i_node = r->i_nodes[i]; k_node = r->k_nodes[i];

        double dcu_p = 0; double dcu_m = 0; 

        if(k_node != NULL){
            dcu_p = fmax(r->v[i], 0.0) * r->prob + fmin(r->v[i], 0.0) * k_node->prob;
        }else{
            dcu_p = fmax(r->v[i], 0.0) * r->prob;
        }
        if(i_node != NULL){
            dcu_m = fmax(i_node->v[i], 0.0)*i_node->prob + fmin(i_node->v[i], 0.0)*r->prob;
        }

        r->dcu -= (G.dt/G.dx[i])*(dcu_p-dcu_m);
    }
}

void update_ctu(TreeNode* r, Grid G){
    if (r == NULL){
        return;
    }
    update_ctu(r->left, G);
    update_ctu(r->right, G);

    for(int i = 0; i < G.dim; i++){
        TreeNode* i_node = r->i_nodes[i];
        TreeNode* j_node; TreeNode* p_node;
        if(i_node!=NULL){
            double F = G.dt*(r->prob-i_node->prob)/(2*G.dx[i]);
            for(int j = 0; j < G.dim; j++){
                if (j!=i){
                    j_node = r->i_nodes[j];
                    p_node = i_node->i_nodes[j];

                    r->ctu[j]      -= fmax(i_node->v[i], 0.0) * fmax(r->v[j], 0.0) * F;
                    i_node->ctu[j] -= fmin(i_node->v[i], 0.0) * fmax(i_node->v[j], 0.0) * F;

                    if(j_node!=NULL){
                        j_node->ctu[j] -= fmax(i_node->v[i], 0.0) * fmin(j_node->v[j], 0.0) * F;
                    }
                    if(p_node!=NULL){
                        p_node->ctu[j] -= fmin(i_node->v[i], 0.0) * fmin(p_node->v[j], 0.0) * F;
                    }
                }
            }

            // High-Resolution Correction Terms
            double th;
            if (i_node->v[i]>0){
                TreeNode* i_i_node = i_node->i_nodes[i];
                if(i_i_node != NULL){
                    th = (i_node->prob-i_i_node->prob)/(r->prob-i_node->prob);
                }else{
                    th = (i_node->prob)/(r->prob-i_node->prob); 
                }
            }else{
                TreeNode* k_node = r->k_nodes[i]; 
                if(k_node != NULL){
                    th = (k_node->prob-r->prob)/(r->prob-i_node->prob);
                }else{
                    th = (-r->prob)/(r->prob-i_node->prob);
                }
            }

            i_node->ctu[i] += fabs(i_node->v[i])*(G.dx[i]/G.dt - fabs(i_node->v[i]))*F*mc(th);
        }
    }
}

void godunov_method(TreeNode** P, Grid G){
    update_dcu(*P, G);
    update_ctu(*P, G);
}

void update_prob(TreeNode* r, Grid G){
    if (r == NULL){
        return;
    }

    update_prob(r->left, G);
    update_prob(r->right, G);

    r->prob += r->dcu;
    for(int i = 0; i < G.dim; i++){
        TreeNode* i_node = r->i_nodes[i]; 
        if(i_node != NULL){
            r->prob -= (G.dt/G.dx[i])*(r->ctu[i]-i_node->ctu[i]);
        }else{
            r->prob -= (G.dt/G.dx[i])*(r->ctu[i]);
        }
    }
    r->prob = fmax(r->prob, 0.0); 
}

void mark_cells(TreeNode* r, Grid G, double* del_probs, uint64_t* del_keys, int* idx){
    if (r == NULL){
        return;
    }
    mark_cells(r->left, G, del_probs, del_keys, idx);
    mark_cells(r->right, G, del_probs, del_keys, idx);

    r->ik_f = 0; bool DELETE = true;
    if (r->prob < G.thresh){

        for(int i = 0; i < G.dim; i++){
            // Looking at Backwards Node
            TreeNode* i_node = r->i_nodes[i]; 

            if(i_node != NULL){
                if ((i_node->v[i]>0)&&(i_node->prob >= G.thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int j = 0; j < G.dim; j++){
                        if(j!=i){
                            TreeNode* i_i_node = i_node->i_nodes[j]; 
                            if(i_i_node != NULL){
                                if ((i_i_node->v[i]>0)&&(i_i_node->v[j]>0)&&(i_i_node->prob >= G.thresh)){
                                    DELETE = false;
                                    break;
                                }
                            }
                            TreeNode* i_k_node = i_node->k_nodes[j]; 
                            if(i_k_node != NULL){
                                if ((i_k_node->v[i]>0)&&(i_k_node->v[j]<0)&&(i_k_node->prob >= G.thresh)){
                                    DELETE = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            // Looking at Forwards Node
            TreeNode* k_node = r->k_nodes[i]; 

            if(k_node != NULL){
                if ((k_node->v[i]<0)&&(k_node->prob >= G.thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int j = 0; j < G.dim; j++){
                        if(j!=i){
                            TreeNode* k_i_node = k_node->i_nodes[j];
                            if(k_i_node != NULL){
                                if ((k_i_node->v[i]<0)&&(k_i_node->v[j]>0)&&(k_i_node->prob >= G.thresh)){
                                    DELETE = false;
                                    break;
                                }
                            }
                            
                            TreeNode* k_k_node = k_node->k_nodes[j];
                            if(k_k_node != NULL){
                                if ((k_k_node->v[i]<0)&&(k_k_node->v[j]<0)&&(k_k_node->prob >= G.thresh)){
                                    DELETE = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        if(DELETE){
            del_probs[*idx] = r->prob; 
            del_keys[*idx] = r->key; 
            *idx += 1; 
        }
    }
}

// REF- different definition of qsort_r in MAC and Linux
#ifdef __linux__ 
  int compare_indices(const void *a, const void *b, void *del_probs) {
#else  
  int compare_indices(void *del_probs, const void *a, const void *b) {
#endif  
    const double *double_list = (double *)del_probs;
    int idx1 = *(const int *)a;
    int idx2 = *(const int *)b;

    if (double_list[idx1] < double_list[idx2])
        return -1;
    else if (double_list[idx1] > double_list[idx2])
        return 1;
    else
        return 0;
}

void sort_by_double(double *del_probs, uint64_t *del_keys, size_t n) {
    int *indices = malloc(n * sizeof(int));
    if (indices == NULL) {
        perror("Error: failed to allocate memory for indices");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < n; i++) {
        indices[i] = i;
    }

    // Sort the indices based on the corresponding doubles    
    // REF- different definition of qsort_r in MAC and Linux
    #ifdef __linux__ 
        qsort_r(indices, n, sizeof(int), compare_indices, del_probs);
    #else  
        qsort_r(indices, n, sizeof(int), del_probs, compare_indices);  
    #endif 

    // Create a temporary array to hold the sorted uint64_t's
    double *sorted_del_probs = malloc(n * sizeof(uint64_t));
    uint64_t *sorted_del_keys = malloc(n * sizeof(uint64_t));
    if ((sorted_del_probs == NULL)||(sorted_del_keys == NULL)) {
        perror("Error: failed to allocate memory for sorted probs or keys");
        free(indices);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < n; i++) {
        sorted_del_probs[i] = del_probs[indices[i]];
        sorted_del_keys[i] = del_keys[indices[i]];
    }

    for (size_t i = 0; i < n; i++) {
        del_probs[i] = sorted_del_probs[i];
        del_keys[i] = sorted_del_keys[i];
    }

    free(indices);
    free(sorted_del_probs);
    free(sorted_del_keys);
}

void delete_cells(TreeNode** P, Grid G, double* del_probs, uint64_t* del_keys, int idx){
    double del_prob_sum = 0; 
    for(int i = 0; i < idx; i++){
        del_prob_sum += del_probs[i]; 
        if(del_probs[i]/(1-del_prob_sum) < G.thresh){
            *P = delete_node(*P, del_keys[i], G);
        }else{
            break;
        }
    }
}

void prune_tree(TreeNode** P, Grid G, int tot_count){
    double* del_probs = malloc(tot_count * sizeof(double));
    uint64_t* del_keys = malloc(tot_count * sizeof(uint64_t));    
    if ((del_keys == NULL)||(del_keys == NULL)){
        perror("Error: memory allocation failure during del_probs/del_keys creation");
    }
    int idx = 0; 
    mark_cells(*P, G, del_probs, del_keys, &idx);
    sort_by_double(del_probs, del_keys, idx);
    delete_cells(P, G, del_probs, del_keys, idx);
    initialize_ik_nodes(*P, *P, &G);
    free(del_probs); 
    free(del_keys); 
}

void meas_up_recursive(void (*h)(double*, double*, double*, double*), TreeNode* r, Grid G, Meas M, Traj T){
    if (r == NULL){
        return;
    }

    meas_up_recursive(h, r->left, G, M, T);
    meas_up_recursive(h, r->right, G, M, T);

    double x[G.dim];
    for(int i = 0; i < G.dim; i++){
        x[i] = G.dx[i]*r->state[i]+G.center[i];
    }

    double y[M.dim];
    (*h)(y, x, G.dx, T.coef); 

    double prob = gauss_probability(M.dim, y, M);   
    r->prob *= prob;
}

void run_gbees(void (*f)(double*, double*, double*, double*), void (*h)(double*, double*, double*, double*), double (*BOUND_f)(double*, double*), Grid G, Meas M, Traj T, char* P_DIR, char* M_DIR, int NUM_DIST, int NUM_MEAS, int DEL_STEP, int OUTPUT_FREQ, int DIM_h, bool OUTPUT, bool RECORD, bool MEASURE, bool BOUNDS){
    char* P_PATH;
    uint64_t max_key; 
    int a_count; 
    int tot_count; 
    double RECORD_TIME = M.T/(NUM_DIST-1);      

    TreeNode* P = NULL; 

    printf("Initializing distribution...\n\n");

    initialize_grid(f, &P, &G, M, T, BOUNDS, BOUND_f); 
    if(BOUNDS){G.lo_bound = DBL_MAX; G.hi_bound = -DBL_MAX; set_bounds(P, &G);} 
    normalize_tree(P, &G, &max_key, &a_count, &tot_count); 

    printf("Entering time marching...\n\n");

    clock_t start = clock(); 
    double tt = 0;
    for(int nm = 0; nm < NUM_MEAS; nm++){
        printf("Timestep: %d-0, Program time: %f s, Sim. time: %f", nm, ((double)(clock()-start))/CLOCKS_PER_SEC, tt); 
        printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", a_count, tot_count, (double)(max_key)/(pow(2,64)-1)*100); 
        if(RECORD){P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", 0); record_data(P, P_PATH, G, tt); free(P_PATH);};

        double mt = 0; int record_count = 1; int step_count = 1; double rt; 
        while(fabs(mt - M.T) > TOL) { // time between discrete measurements

            rt = 0;
            while (rt < RECORD_TIME) { // time between PDF recordings

                grow_tree(f, &P, G, T, BOUNDS, BOUND_f);
                check_cfl_condition(P, &G); 
                G.dt = fmin(G.dt, RECORD_TIME - rt);
                rt += G.dt;
                godunov_method(&P, G);
                update_prob(P, G);
                normalize_tree(P, &G, &max_key, &a_count, &tot_count); 

                if (step_count % DEL_STEP == 0) { // deletion procedure
                    prune_tree(&P, G, tot_count);
                    normalize_tree(P, &G, &max_key, &a_count, &tot_count); 
                }

                if ((OUTPUT) && (step_count % OUTPUT_FREQ == 0)) { // print size to terminal
                    printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt); 
                    printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", a_count, tot_count, (double)(max_key)/(pow(2,64)-1)*100); 
                }

                step_count += 1; G.dt = DBL_MAX;
            }
            
            if (((step_count-1) % OUTPUT_FREQ != 0)||(!OUTPUT)){ // print size to terminal  
                printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count - 1, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt); 
                printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", a_count, tot_count, (double)(max_key)/(pow(2,64)-1)*100); 
            }
            
            if (RECORD) { // record PDF
                P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", record_count); 
                record_data(P, P_PATH, G, tt + mt + rt);
                record_count += 1;
                free(P_PATH);
            }

            mt += rt;
        }

        tt += mt;
        if ((MEASURE) && (nm < NUM_MEAS - 1)) { // peform discrete measurement update

            printf("\nPERFORMING BAYESIAN UPDATE AT: %f TU...\n\n", tt);

            // freeing previous measurement 
            Meas_free(&M);                                                    

            // reading new measurement
            char* M_FILE_NAME = "measurement";                                             
            char* M_FILE_EXT = ".txt"; 
            int length = snprintf(NULL, 0, "%s%d%s", M_FILE_NAME, nm+1, M_FILE_EXT) + 1;
            char* M_FILE = (char *)malloc(length);
            if (M_FILE == NULL) {
                perror("Error: memory allocation failure during discrete measurement update");
                exit(EXIT_FAILURE);
            }
            snprintf(M_FILE, length, "%s%d%s", M_FILE_NAME, nm+1, M_FILE_EXT);
            M = Meas_create(DIM_h, M_DIR, M_FILE);                                       
            free(M_FILE); 

            // performing discrete update
            meas_up_recursive(h, P, G, M, T);                                     
            normalize_tree(P, &G, &max_key, &a_count, &tot_count); 
            prune_tree(&P, G, tot_count);
            normalize_tree(P, &G, &max_key, &a_count, &tot_count); 
        }
    }
    
    Meas_free(&M); 
    Grid_free(&G); 
    Traj_free(&T);
    Tree_free(P);

    printf("Time marching complete.\n");

    return;
}
