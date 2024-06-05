/*==============================================================================

MAIN.CPP

==============================================================================*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <string>
#include <armadillo>

using namespace arma;

const int DIM = 6;

class Traj{
public:
    double mu;
    double T;
};

void Record_Data(const std::string& filename, const vec& mean, const mat& cov, double t){
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

vec f_classical(const vec& mean_1, double dt, Traj t){
    vec f1, f2, f3, f4, mean_2, mean_3, mean_4;
    f1 = {0, 0, 0, 0, 0, pow(t.mu/pow(mean_1(0),3),0.5)};
    mean_2 = mean_1 + 0.5*dt*f1;
    f2 = {0, 0, 0, 0, 0, pow(t.mu/pow(mean_2(0),3),0.5)};
    mean_3 = mean_1 + 0.5*dt*f2;
    f3 = {0, 0, 0, 0, 0, pow(t.mu/pow(mean_3(0),3),0.5)};
    mean_4 = mean_1 + dt*f3;
    f4 = {0, 0, 0, 0, 0, pow(t.mu/pow(mean_4(0),3),0.5)};
    return mean_1 + (dt/6.0)*(f1 + 2*f2 + 2*f3 + f4);
}

mat F_classical(const vec& mean, double dt, Traj t) {
    double epsilon = 1e-6;
    vec fx = f_classical(mean, dt, t);

    mat J(DIM, DIM);
    vec mean_eps = mean;

    for (size_t j = 0; j < DIM; ++j) {
        mean_eps = mean;
        mean_eps(j) += epsilon;
        vec fx_eps = f_classical(mean_eps, dt, t);

        J.col(j) = (fx_eps - fx) / epsilon;
    }

    return J;
}

vec f_equinoctal(const vec& mean_1, double dt, Traj t){
    vec f1, f2, f3, f4, mean_2, mean_3, mean_4;
    f1 = {0, 0, 0, 0, 0, pow(t.mu*mean_1(0),0.5)*pow((1+mean_1(1)*cos(mean_1(5))+mean_1(2)*sin(mean_1(5)))/mean_1(0),2.0)};
    mean_2 = mean_1 + 0.5*dt*f1;
    f2 = {0, 0, 0, 0, 0, pow(t.mu*mean_2(0),0.5)*pow((1+mean_2(1)*cos(mean_2(5))+mean_2(2)*sin(mean_2(5)))/mean_2(0),2.0)};
    mean_3 = mean_1 + 0.5*dt*f2;
    f3 = {0, 0, 0, 0, 0, pow(t.mu*mean_3(0),0.5)*pow((1+mean_3(1)*cos(mean_3(5))+mean_3(2)*sin(mean_3(5)))/mean_3(0),2.0)};
    mean_4 = mean_1 + dt*f3;
    f4 = {0, 0, 0, 0, 0, pow(t.mu*mean_4(0),0.5)*pow((1+mean_4(1)*cos(mean_4(5))+mean_4(2)*sin(mean_4(5)))/mean_4(0),2.0)};
    return mean_1 + (dt/6.0)*(f1 + 2*f2 + 2*f3 + f4);
}

mat F_equinoctal(const vec& mean, double dt, Traj t) {
    double epsilon = 1e-6;
    vec fx = f_equinoctal(mean, dt, t);

    mat J(DIM, DIM);
    vec mean_eps = mean;

    for (size_t j = 0; j < DIM; ++j) {
        mean_eps = mean;
        mean_eps(j) += epsilon;
        vec fx_eps = f_equinoctal(mean_eps, dt, t);

        J.col(j) = (fx_eps - fx) / epsilon;
    }

    return J;
}

vec f_cartesian(const vec& mean_1, double dt, Traj t){
    vec f1, f2, f3, f4, mean_2, mean_3, mean_4;
    f1 = {mean_1(3), mean_1(4), mean_1(5), -(t.mu/pow(pow(mean_1(0),2)+pow(mean_1(1),2)+pow(mean_1(2),2),1.5))*mean_1(0), -(t.mu/pow(pow(mean_1(0),2)+pow(mean_1(1),2)+pow(mean_1(2),2),1.5))*mean_1(1), -(t.mu/pow(pow(mean_1(0),2)+pow(mean_1(1),2)+pow(mean_1(2),2),1.5))*mean_1(2)};
    mean_2 = mean_1 + 0.5*dt*f1;
    f2 = {mean_2(3), mean_2(4), mean_2(5), -(t.mu/pow(pow(mean_2(0),2)+pow(mean_2(1),2)+pow(mean_2(2),2),1.5))*mean_2(0), -(t.mu/pow(pow(mean_2(0),2)+pow(mean_2(1),2)+pow(mean_2(2),2),1.5))*mean_2(1), -(t.mu/pow(pow(mean_2(0),2)+pow(mean_2(1),2)+pow(mean_2(2),2),1.5))*mean_2(2)};
    mean_3 = mean_1 + 0.5*dt*f2;
    f3 = {mean_3(3), mean_3(4), mean_3(5), -(t.mu/pow(pow(mean_3(0),2)+pow(mean_3(1),2)+pow(mean_3(2),2),1.5))*mean_3(0), -(t.mu/pow(pow(mean_3(0),2)+pow(mean_3(1),2)+pow(mean_3(2),2),1.5))*mean_3(1), -(t.mu/pow(pow(mean_3(0),2)+pow(mean_3(1),2)+pow(mean_3(2),2),1.5))*mean_3(2)};
    mean_4 = mean_1 + dt*f3;
    f4 = {mean_4(3), mean_4(4), mean_4(5), -(t.mu/pow(pow(mean_4(0),2)+pow(mean_4(1),2)+pow(mean_4(2),2),1.5))*mean_4(0), -(t.mu/pow(pow(mean_4(0),2)+pow(mean_4(1),2)+pow(mean_4(2),2),1.5))*mean_4(1), -(t.mu/pow(pow(mean_4(0),2)+pow(mean_4(1),2)+pow(mean_4(2),2),1.5))*mean_4(2)};
    return mean_1 + (dt/6.0)*(f1 + 2*f2 + 2*f3 + f4);
}

mat F_cartesian(const vec& mean, double dt, Traj t) {
    double epsilon = 1e-6;
    vec fx = f_cartesian(mean, dt, t);

    mat J(DIM, DIM);
    vec mean_eps = mean;

    for (size_t j = 0; j < DIM; ++j) {
        mean_eps = mean;
        mean_eps(j) += epsilon;
        vec fx_eps = f_cartesian(mean_eps, dt, t);

        J.col(j) = (fx_eps - fx) / epsilon;
    }

    return J;
}

int main(){
    //===================================== Begin User Input ====================================================
    std::cout << "Reading in user inputs...\n" << std::endl;

    bool OUTPUT = true;     // Print info to terminal
    bool RECORD = true;     // Write EKF to .txt file (for epochs)
    int num_dist = 4;       // Number of distributions recorded per revolution
    const int REV = 1;      // Number of revolutions
    int oetype = 2;         // 0: classical orbit elements, 1: equinoctal orbit elements, 2: cartesian position and velocity
    int body = 0;           // 0: Earth, 1: Europa

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
    vec mean(DIM, fill::zeros);
    mat cov(DIM, DIM, fill::zeros);
    mat Q(DIM, DIM, fill::zeros);

    std::cout << "Reading in discrete measurements..." << std::endl; std::cout << std::endl;

    std::ifstream measurement_file(FILE_PATH + "/measurements.txt"); // Measurement file
    std::string line;
    int r = 0;
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < 1)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            mean(r,c) = std::stod(line);
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
            cov(r,c) = std::stod(line);
            c++;
        }
        r++;
    }

    Traj lyap{}; // Trajectory object
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.mu = std::stod(line); // Read in mu
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.T = REV*std::stod(line);  // Read in T

    vec tspan = linspace<vec>(0, REV*lyap.T, (REV*num_dist)+1);

    measurement_file.close();
    //===========================================================================================================
    std::cout << "Initializing Distribution..." << std::endl; std::cout << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed{};

    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
    std::cout << "Program time: " << elapsed.count() << " s, Sim. time: " << 0 << " s" << std::endl;
    if(RECORD){std::string filename = FILE_PATH + "/Revs" + std::to_string(REV) + "/pdf_0.txt"; Record_Data(filename, mean, cov, 0);}

    double h = 0.1; int record_count = 1;
    for(int i=1; i < tspan.n_elem; i++){
        
        if(oetype==0){
            mean = f_classical(mean, h, lyap);
            cov  = F_classical(mean, h, lyap) * cov * trans(F_classical(mean, h, lyap)) + Q;
        }else if(oetype==1){
            mean = f_equinoctal(mean, h, lyap);
            cov  = F_equinoctal(mean, h, lyap) * cov * trans(F_equinoctal(mean, h, lyap)) + Q;
        }else if(oetype==2){
            mean = f_cartesian(mean, h, lyap);
            cov  = F_cartesian(mean, h, lyap) * cov * trans(F_equinoctal(mean, h, lyap)) + Q;
        }

        if(OUTPUT){
            finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
            std::cout << "Program time: " << elapsed.count() << " s, Sim. time: " << tspan(i) << " s" << std::endl;
        }

        if(RECORD){ // Record PF Scatter Plot
            std::string filename = FILE_PATH + "/Revs" + std::to_string(REV) + "/pdf_" + std::to_string(record_count) + ".txt"; Record_Data(filename, mean, cov, tspan(i));
            record_count+=1;
        }
    }
    return 0;
}