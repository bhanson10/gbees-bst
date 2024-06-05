/*==============================================================================

MAIN.CPP

==============================================================================*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <armadillo>

using namespace arma;

const int DIM = 6;

// Runge-Kutta 8(7) coefficients
const vec c1 = {1.0/18, 1.0/12, 1.0/8, 5.0/16, 3.0/8, 59.0/400, 93.0/200, 5490023248.0/9719169821, 13.0/20, 1201146811.0/1299019798, 1.0, 1.0};
mat aT = {{1.0/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {1.0/48, 1.0/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {1.0/32, 0, 3.0/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {5.0/16, 0, -75.0/64, 75.0/64, 0, 0, 0, 0, 0, 0, 0, 0, 0},
          {3.0/80, 0, 0, 3.0/16, 3.0/20, 0, 0, 0, 0, 0, 0, 0, 0},
          {29443841.0/614563906, 0, 0, 77736538.0/692538347, -28693883.0/1125000000, 23124283.0/1800000000, 0, 0, 0, 0, 0, 0, 0},
          {16016141.0/946692911, 0, 0, 61564180.0/158732637, 22789713.0/633445777, 545815736.0/2771057229, -180193667.0/1043307555, 0, 0, 0, 0, 0, 0},
          {39632708.0/573591083, 0, 0, -433636366.0/683701615, -421739975.0/2616292301, 100302831.0/723423059, 790204164.0/839813087, 800635310.0/3783071287, 0, 0, 0, 0, 0},
          {246121993.0/1340847787, 0, 0, -37695042795.0/15268766246, -309121744.0/1061227803, -12992083.0/490766935, 6005943493.0/2108947869, 393006217.0/1396673457, 123872331.0/1001029789, 0, 0, 0, 0},
          {-1028468189.0/846180014, 0, 0, 8478235783.0/508512852, 1311729495.0/1432422823, -10304129995.0/1701304382, -48777925059.0/3047939560, 15336726248.0/1032824649, -45442868181.0/3398467696, 3065993473.0/597172653, 0, 0, 0},
          {185892177.0/718116043, 0, 0, -3185094517.0/667107341, -477755414.0/1098053517, -703635378.0/230739211, 5731566787.0/1027545527, 5232866602.0/850066563, -4093664535.0/808688257, 3962137247.0/1805957418, 65686358.0/487910083, 0, 0},
          {403863854.0/491063109, 0, 0, -5068492393.0/434740067, -411421997.0/543043805, 652783627.0/914296604, 11173962825.0/925320556, -13158990841.0/6184727034, 3936647629.0/1978049680, -160528059.0/685178525, 248638103.0/1413531060, 0, 0}};
const mat a = trans(aT);
const vec b7 = {14005451.0/335480064, 0.0, 0.0, 0.0, 0.0, -59238493.0/1068277825, 181606767.0/758867731, 561292985.0/797845732, -1041891430.0/1371343529, 760417239.0/1151165299, 118820643.0/751138087, -528747749.0/2220607170, 1.0/4};
const vec b8 = {13451932.0/455176623, 0.0, 0.0, 0.0, 0.0, -808719846.0/976000145, 1757004468.0/5645159321, 656045339.0/265891186, -3867574721.0/1518517206, 465885868.0/322736535, 53011238.0/667516719, 2.0/45, 0.0};
const double EPS = 2.220446049250313e-16;
const double POW = 1.0/8;

struct Member {
    vec state;
};

class Traj{
public:
    double mu;
    double T;
};

vec classical(double t, vec y, double mu){
    return {0, 0, 0, 0, 0, sqrt(mu/pow(y(0),3.0))};
}

vec equinoctal(double t, vec y, double mu){
    return {0, 0, 0, 0, 0, pow(mu*y(0),0.5)*pow((1+y(1)*cos(y(5))+y(2)*sin(y(5)))/y(0),2.0)};
}

vec cartesian(double t, vec y, double mu){
    return {y(3), y(4), y(5), -(mu/pow(pow(y(0),2)+pow(y(1),2)+pow(y(2),2),1.5))*y(0), -(mu/pow(pow(y(0),2)+pow(y(1),2)+pow(y(2),2),1.5))*y(1), -(mu/pow(pow(y(0),2)+pow(y(1),2)+pow(y(2),2),1.5))*y(2)};
};

class Ensemble{
public:
    vec mean;
    mat cov;
    int size;
    std::vector<Member> members;

    explicit Ensemble(int n, int s) : mean(n, fill::zeros), cov(n, n, fill::zeros), size(s) {}

    void Initialize_Ensemble(){
        members.clear();

        mat samples = mvnrnd(mean, cov, size);

        for(int q = 0; q < size; q++){
            Member m{};
            m.state = samples.col(q);

            members.push_back(m);
        }
    }

    void Record_Data(const std::string& filename, double t){
        std::ofstream myfile; myfile.open(filename);
        myfile << t << std::endl << std::endl;

        for(int i=0; i < DIM; i++){
            myfile << mean(i) << " ";
        }
        myfile << std::endl << std::endl;

        for(int i=0; i < DIM; i++){
            for(int j=0; j < DIM; j++){
                myfile << cov(i,j) << " ";
            }
            myfile << std::endl;
        }
    }

    void RK87(double t0, double t1, double h, Traj lyap, int oetype){
        for(Member& m : members) {
            vec y = m.state;

            double tol = 1e-6;
            double hmax = 5;
            int n_reject = 0;
            int reject = 0;
            mat f(DIM,13,fill::zeros);
            double t = t0;
            double tfinal = t1;
            while (t < tfinal) {
                if ((t + h) > tfinal) {
                    h = tfinal - t;
                }

                if (oetype == 0) {
                    f.col(0) = classical(t, y, lyap.mu);
                } else if (oetype == 1) {
                    f.col(0) = equinoctal(t, y, lyap.mu);
                } else if (oetype == 2) {
                    f.col(0) = cartesian(t, y, lyap.mu);
                }
                for (int j = 0; j <= 11; j++) {
                    vec yh = y + h * f * a.col(j);
                    if (oetype == 0) {
                        f.col(j + 1) = classical(t + c1(j) * h, yh, lyap.mu);
                    } else if (oetype == 1) {
                        f.col(j + 1) = equinoctal(t + c1(j) * h, yh, lyap.mu);
                    } else if (oetype == 2) {
                        f.col(j + 1) = cartesian(t + c1(j) * h, yh, lyap.mu);
                    }
                }

                vec sol2 = y + h * f * b8;
                vec sol1 = y + h * f * b7;

                double error_1 = norm(sol1 - sol2);
                double error_step = abs(error_1);
                double tau = tol * std::max(norm(y, "inf"), 1.0);

                if (error_step <= tau) {
                    t = t + h;
                    y = sol2;
                    reject--;
                } else {
                    n_reject++;
                    reject = 1;
                }

                if (error_step == 0.0) {
                    error_step = EPS * 10.0;
                }
                h = std::min(hmax, 0.9 * h * pow((tau / error_step), POW));
                if (abs(h) <= EPS) {
                    if (reject == 0) {
                        std::cout << "WARNING!!! ode87. Step is very small!!!" << std::endl;
                        h = EPS * 100;
                    } else {
                        std::cout << "Error in ode87. Step is too small." << std::endl;
                    }
                }
            }
            m.state = y;
        }
    }

    vec perturb(const mat& Q){
        vec sum(DIM,fill::zeros);
        for(auto & m : members) {
            m.state = mvnrnd(m.state, Q, 1);
            sum += m.state;
        }
        return sum;
    }
};

int main(){
    //===================================== Begin User Input ====================================================
    std::cout << "Reading in user inputs...\n" << std::endl;

    bool OUTPUT = true; // Print info to terminal
    bool RECORD = true; // Write Ensemble to .txt file (for epochs)
    int num_dist = 4;  // Number of distributions recorded per revolution
    const int REV = 6;  // Number of revolutions
    int oetype = 2;     // 0: classical orbit elements, 1: equinoctal orbit elements, 2: cartesian position and velocity
    int body = 1;       // 0: Earth, 1: Europa

    std::string FILE_PATH;
    if(oetype==0) {
        if(body==0){
            FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Earth/Classical";
        }else if(body==1) {
            FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Classical";
        }
    }
    else if(oetype==1){
        if(body==0){
            FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Earth/Equinoctial";
        }else if(body==1) {
            FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Equinoctial";
        }
    }
    else if(oetype==2){
        if(body==0){
            FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Earth/Cartesian";
        }else if(body==1) {
            FILE_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Cartesian";
        }
    }
    //===========================================================================================================

    //===================================== Read in measurement/trajectory info =================================
    const int SIZE = 10000;
    mat Q(DIM, DIM, fill::zeros);
    Ensemble M(DIM, SIZE);

    std::cout << "Reading in discrete measurements..." << std::endl; std::cout << std::endl;

    std::ifstream measurement_file(FILE_PATH + "/measurements.txt"); // Measurement file
    std::string line;
    int r = 0;
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < 1)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            M.mean(r,c) = std::stod(line);
            r++;
        }
        c++;
    }

    r = 0;
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < DIM)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            M.cov(r,c) = std::stod(line);
            c++;
        }
        r++;
    }

    Traj lyap{}; // Trajectory object
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.mu = std::stod(line); // Read in mu
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.T = std::stod(line);  // Read in T

    vec tspan = linspace<vec>(0, REV*lyap.T, (REV*num_dist)+1);

    measurement_file.close();
    //===========================================================================================================
    std::cout << "Initializing Distribution..." << std::endl; std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed{};

    M.Initialize_Ensemble();

    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
    std::cout << "Program time: " << elapsed.count() << " s, Sim. time: " << 0 << " s, Members: " << M.size << std::endl;
    if(RECORD){std::string filename = FILE_PATH + "/Revs" + std::to_string(REV) + "/pdf_0.txt"; M.Record_Data(filename, 0);}

    double h = 0.1; int record_count = 1;
    for(int i=1; i < tspan.n_elem; i++){
        M.RK87(tspan(i-1), tspan(i), h, lyap, oetype);

        M.mean = M.perturb(Q)/SIZE;

        int ens_count = 0;
        mat delX(DIM, SIZE, fill::zeros);
        for (Member& m : M.members){
            vec del_xj = m.state - M.mean;
            delX.col(ens_count) = del_xj;
            ens_count++;
        }
        M.cov = delX*trans(delX)/(SIZE-1);

        if(OUTPUT){
            finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
            std::cout << "Program time: " << elapsed.count() << " s, Sim. time: " << tspan(i) << " s, Members: " << M.size << std::endl;
        }

        if(RECORD){ // Record PF Scatter Plot
            std::string filename = FILE_PATH + "/Revs" + std::to_string(REV) + "/pdf_" + std::to_string(record_count) + ".txt"; M.Record_Data(filename, tspan(i));
            record_count+=1;
        }
    }
    return 0;
}