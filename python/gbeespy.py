import ctypes as ct
import time 
import sys
sys.setrecursionlimit(10000)

lib = ct.CDLL("./gbees.so")
TOL = 1E-8

class Meas(ct.Structure):
    _fields_ = [
        ("mean", ct.POINTER(ct.c_double)),
        ("cov", ct.POINTER(ct.POINTER(ct.c_double))),
        ("T", ct.c_double)
    ]

def Meas_create(dim, M_DIR, M_FILE):
    lib.Meas_create.argtypes = [ct.c_int, ct.c_char_p, ct.c_char_p]
    lib.Meas_create.restype = Meas

    M_DIR = ct.c_char_p(M_DIR.encode('utf-8'))
    M_FILE = ct.c_char_p(M_FILE.encode('utf-8'))
    return lib.Meas_create(dim, M_DIR, M_FILE)

class Grid(ct.Structure):
    _fields_ = [
        ("dim", ct.c_int),
        ("thresh", ct.c_double),
        ("dt", ct.c_double),
        ("center", ct.POINTER(ct.c_double)),
        ("dx", ct.POINTER(ct.c_double))
    ]

def Grid_create(dim, thresh, center, dx):
    lib.Grid_create.argtypes = [ct.c_int, ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
    lib.Grid_create.restype = Grid

    dx = (ct.c_double * dim)(*dx)
    return lib.Grid_create(dim, thresh, center, dx)

class Traj(ct.Structure):
    _fields_ = [
        ("coef", ct.POINTER(ct.c_double))
    ]

def Traj_create(n, coef):
    lib.Traj_create.argtypes = [ct.c_int, ct.POINTER(ct.c_double)]
    lib.Traj_create.restype = Traj

    coef = (ct.c_double * n)(*coef)
    return lib.Traj_create(n, coef)

class TreeNode(ct.Structure):
    pass

TreeNode._fields_ = [
    ("key", ct.c_ulong),
    ("prob", ct.c_double),
    ("v", ct.POINTER(ct.c_double)),
    ("u", ct.POINTER(ct.c_double)),
    ("w", ct.POINTER(ct.c_double)),
    ("ctu", ct.POINTER(ct.c_double)),
    ("state", ct.POINTER(ct.c_int)),
    ("i_nodes", ct.POINTER(ct.POINTER(TreeNode))),
    ("k_nodes", ct.POINTER(ct.POINTER(TreeNode))),
    ("dcu", ct.c_double),
    ("cfl_dt", ct.c_double),
    ("new_f", ct.c_int),
    ("ik_f", ct.c_int),
    ("del_f", ct.c_int),
    ("left", ct.POINTER(TreeNode)),
    ("right", ct.POINTER(TreeNode))
]

def initialize_vuw(r, G, T, ADV):
    if not r:
        return

    initialize_vuw(r.contents.left, G, T, ADV)
    initialize_vuw(r.contents.right, G, T, ADV)

    if r.contents.new_f == 0:
        x = [None] * G.dim
        for i in range(G.dim):
            x[i] = G.dx[i] * r.contents.state[i] + G.center[i]

        v = ADV(x, G.dx, T)

        for i in range(G.dim):
            r.contents.v[i] = v[i]
            r.contents.u[i] = min(v[i], 0.0)
            r.contents.w[i] = max(v[i], 0.0)
        r.contents.new_f = 1

        sum = 0
        for q in range(G.dim):
            sum += abs(r.contents.v[q]) / G.dx[q]
        r.contents.cfl_dt = 1 / sum

def run_gbees(ADV, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE):
    lib.initialize_grid.argtypes = [ct.POINTER(Grid), Meas, Traj]
    lib.initialize_grid.restype = ct.POINTER(TreeNode)
    lib.normalize_tree.argtypes = [ct.POINTER(TreeNode)]
    lib.normalize_tree.restype = None
    lib.get_tree_info.argtypes = [ct.POINTER(TreeNode), ct.POINTER(Grid), ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
    lib.get_tree_info.restype = None
    lib.record_data.argtypes = [ct.POINTER(TreeNode), ct.c_char_p, Grid, ct.c_double]
    lib.record_data.restype = None
    lib.grow_tree.argtypes = [ct.POINTER(ct.POINTER(TreeNode)), Grid]
    lib.grow_tree.restype = None
    lib.update_prob.argtypes = [ct.POINTER(TreeNode), ct.POINTER(Grid), ct.c_double]
    lib.update_prob.restype = None
    lib.prune_tree.argtypes = [ct.POINTER(ct.POINTER(TreeNode)), Grid]
    lib.prune_tree.restype = None
    lib.measurement_update_recursive.argtypes = [ct.POINTER(TreeNode), Grid, Meas]
    lib.measurement_update_recursive.restype = None
    lib.Meas_free.argtypes = [ct.POINTER(Meas), ct.c_int]
    lib.Meas_free.restype =  None
    lib.Grid_free.argtypes = [ct.POINTER(Grid)]
    lib.Grid_free.restype = None
    lib.Traj_free.argtypes = [ct.POINTER(Traj)]
    lib.Traj_free.restype = None
    lib.Tree_free.argtypes = [ct.POINTER(TreeNode)]
    lib.Tree_free.restype = None

    RECORD_TIME = M.T/(NUM_DIST-1)

    print("Initializing distribution...\n\n")
    
    P = lib.initialize_grid(ct.pointer(G), M, T)
    initialize_vuw(P, G, T, ADV)
    lib.normalize_tree(P)

    print("Entering time marching...\n\n")

    start = time.time()
    tt = 0.0
    for nm in range(NUM_MEAS):
        max_key = ct.c_uint64(0)
        a_count = ct.c_int(0)
        tot_count = ct.c_int(0)
        lib.get_tree_info(P, ct.pointer(G), ct.pointer(max_key), ct.pointer(a_count), ct.pointer(tot_count))
        finish = time.time()
        if OUTPUT:
            print("Timestep: " + str(nm) + "-0, Program time: " + str(finish - start) + " s, Sim. time: " + str(tt) + " TU, Active/Total Cells: " + str(int(a_count.value)) + "/" + str(int(tot_count.value)) + ", Max key %: " + str(float(max_key.value)/((2**64)-1)*100) + "%")
        if RECORD:
            P_PATH = P_DIR + "/P" + str(nm) + "/pdf_0.txt"
            lib.record_data(P, P_PATH.encode('utf-8'), G, ct.c_double(0.0))

        mt = 0.0
        record_count = 1
        step_count = 1
        while(abs(mt - M.T) > TOL):

            rt = 0.0
            while(rt < RECORD_TIME):
                
                lib.grow_tree(ct.pointer(P), G)
                initialize_vuw(P, G, T, ADV)
                lib.update_prob(P, ct.pointer(G), ct.c_double(RECORD_TIME - rt))

                if (step_count % DEL_STEP == 0):
                    lib.prune_tree(ct.pointer(P), G)

                if OUTPUT and step_count % OUTPUT_FREQ == 0:
                    max_key = ct.c_uint64(0)
                    a_count = ct.c_int(0)
                    tot_count = ct.c_int(0)
                    lib.get_tree_info(P, ct.pointer(G), ct.pointer(max_key), ct.pointer(a_count), ct.pointer(tot_count))
                    finish = time.time()
                    print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(int(a_count.value)) + "/" + str(int(tot_count.value)) + ", Max key %: " + str(float(max_key.value)/((2**64)-1)*100) + "%")
                step_count += 1

                rt += G.dt

            if OUTPUT and step_count % OUTPUT_FREQ != 0:
                max_key = ct.c_uint64(0)
                a_count = ct.c_int(0)
                tot_count = ct.c_int(0)
                lib.get_tree_info(P, ct.pointer(G), ct.pointer(max_key), ct.pointer(a_count), ct.pointer(tot_count))
                finish = time.time()
                print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(int(a_count.value)) + "/" + str(int(tot_count.value)) + ", Max key %: " + str(float(max_key.value)/((2**64)-1)*100) + "%")

            if RECORD:
                if OUTPUT:
                    print("\nRECORDING PDF AT: " + str(tt + mt + rt) + " TU...\n")
                P_PATH = P_DIR + "/P" + str(nm) + "/pdf_" + str(record_count) + ".txt"
                lib.record_data(P, P_PATH.encode('utf-8'), G, ct.c_double(tt + mt + rt))
                record_count += 1
        
            mt += rt

        tt += mt
        if MEASURE and nm < NUM_MEAS-1:
            if OUTPUT:
                print("\nPERFORMING BAYESIAN UPDATE AT: " + str(tt) + " TU...\n")

            M_FILE = "/measurement" + str(nm) + ".txt"
            M = Meas_create(DIM, M_DIR, M_FILE)
            
            # performing discrete update
            lib.measurement_update_recursive(P, G, M);                                     
            lib.normalize_tree(P)
            lib.prune_tree(ct.pointer(P), G)

    lib.Meas_free(ct.pointer(M), G.dim)
    lib.Grid_free(ct.pointer(G))
    lib.Traj_free(ct.pointer(T))
    lib.Tree_free(P)

    print("Time marching complete.")