import ctypes as ct

lib = ct.CDLL("../../gbees.so")
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
        ("dx", ct.POINTER(ct.c_double)),
        ("hi_bound", ct.POINTER(ct.c_double)),
        ("lo_bound", ct.POINTER(ct.c_double))
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
    ("ctu", ct.POINTER(ct.c_double)),
    ("state", ct.POINTER(ct.c_int)),
    ("i_nodes", ct.POINTER(ct.POINTER(TreeNode))),
    ("k_nodes", ct.POINTER(ct.POINTER(TreeNode))),
    ("dcu", ct.c_double),
    ("cfl_dt", ct.c_double),
    ("new_f", ct.c_int),
    ("ik_f", ct.c_int),
    ("del_f", ct.c_int),
    ("bound_val", ct.c_double),
    ("left", ct.POINTER(TreeNode)),
    ("right", ct.POINTER(TreeNode))
]

def run_gbees(ADV_F, BOUND_F, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE, BOUNDS):

    c_P_DIR = P_DIR.encode('utf-8')
    c_P_DIR = ct.create_string_buffer(c_P_DIR)
    c_P_DIR = ct.cast(c_P_DIR, ct.POINTER(ct.c_char))

    c_M_DIR = M_DIR.encode('utf-8')
    c_M_DIR = ct.create_string_buffer(c_M_DIR)
    c_M_DIR = ct.cast(c_M_DIR, ct.POINTER(ct.c_char))

    ADV_CALLBACK_FUNC = ct.CFUNCTYPE(None, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double))
    BOUND_CALLBACK_FUNC = ct.CFUNCTYPE(ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double))

    @ADV_CALLBACK_FUNC
    def c_ADV_F(x, dx, coef):
        v = ADV_F(x, dx, coef)
        for i in range(DIM):
            x[i] = v[i]
    
    if(BOUNDS):
        @BOUND_CALLBACK_FUNC
        def c_BOUND_F(x, coef):
            J = BOUND_F(x, coef)
            return J
        
        lib.run_gbees.argtypes = [ADV_CALLBACK_FUNC, BOUND_CALLBACK_FUNC, Grid, Meas, Traj, ct.c_char_p, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_bool, ct.c_bool, ct.c_bool, ct.c_bool]
        lib.run_gbees.restype = None
    else:
        c_BOUND_F = None
        lib.run_gbees.argtypes = [ADV_CALLBACK_FUNC, ct.c_char_p, Grid, Meas, Traj, ct.c_char_p, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_bool, ct.c_bool, ct.c_bool, ct.c_bool]
        lib.run_gbees.restype = None

    lib.run_gbees(c_ADV_F, c_BOUND_F, G, M, T, c_P_DIR, c_M_DIR, ct.c_int(NUM_DIST), ct.c_int(NUM_MEAS), ct.c_int(DEL_STEP), ct.c_int(OUTPUT_FREQ), ct.c_int(DIM), ct.c_bool(OUTPUT), ct.c_bool(RECORD), ct.c_bool(MEASURE), ct.c_bool(BOUNDS))