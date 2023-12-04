/*==============================================================================

MAIN.CPP

==============================================================================*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <random>
#include <string>
const int DIM = 4; 

struct Particle {
    double x, y, vx, vy;
};

class Traj{
    public: 
        double mu; 
        double T;
};

class MC{
    public: 
        std::array<double,DIM> epoch;
        std::array<double,DIM> std;
        int size; 
        std::vector<Particle> particles; 
        std::vector<std::vector<double>> k1; 
        std::vector<std::vector<double>> k2; 
        std::vector<std::vector<double>> k3; 
        std::vector<std::vector<double>> k4; 

        void Initialize_MC(Traj lyap){
            particles.clear(); 

            std::default_random_engine generator;
            std::normal_distribution<double> distribution1(epoch[0],std[0]);
            std::normal_distribution<double> distribution2(epoch[1],std[1]);
            std::normal_distribution<double> distribution3(epoch[2],std[2]);
            std::normal_distribution<double> distribution4(epoch[3],std[3]);

            for(int q = 0; q < size; q++){
                Particle p;
                p.x = distribution1(generator);
                p.y = distribution2(generator);
                p.vx = distribution3(generator);
                p.vy = distribution4(generator);

                particles.push_back(p); 
            }
        }

        void Record_Data(std::string filename, double t){
            std::ofstream myfile; myfile.open(filename);
            myfile << t << std::endl; 

            for (Particle& p : particles){
                myfile << p.x << " " << p.y << " " << p.vx << " " << p.vy << std::endl; 
            }
        }

        void update_prob(Traj lyap, double dt){
            std::vector<Particle> k1, k2, k3, k4;

            for (Particle& p : particles){ // Calculate k1
                double x = p.x; double y = p.y; double vx = p.vx; double vy = p.vy;

                double r1 = pow(pow(x+lyap.mu,2)+pow(y,2),1.5); 
                double r2 = pow(pow(x-1+lyap.mu,2)+pow(y,2),1.5); 
                double ax = 2*vy+x-(lyap.mu*(x-1+lyap.mu)/r2)-((1-lyap.mu)*(x+lyap.mu)/r1); 
                double ay = -2*vx+y-(lyap.mu*y/r2)-((1-lyap.mu)*y/r1);

                k1.push_back({vx, vy, ax, ay});
            }

            for (int i = 0; i < particles.size(); ++i){ // Calculate k2
                Particle& p = particles[i];
                Particle p_k2 = {p.x + 0.5 * dt * k1[i].x, p.y + 0.5 * dt * k1[i].y, p.vx + 0.5 * dt * k1[i].vx, p.vy + 0.5 * dt * k1[i].vy};
                
                double x = p_k2.x; double y = p_k2.y; double vx = p_k2.vx; double vy = p_k2.vy;

                double r1 = pow(pow(x+lyap.mu,2)+pow(y,2),1.5); 
                double r2 = pow(pow(x-1+lyap.mu,2)+pow(y,2),1.5); 
                double ax = 2*vy+x-(lyap.mu*(x-1+lyap.mu)/r2)-((1-lyap.mu)*(x+lyap.mu)/r1); 
                double ay = -2*vx+y-(lyap.mu*y/r2)-((1-lyap.mu)*y/r1);

                k2.push_back({p.vx + 0.5 * dt * ax, p.vy + 0.5 * dt * ay, ax, ay});
            }

            for (int i = 0; i < particles.size(); ++i){ // Calculate k3
                Particle& p = particles[i];
                Particle p_k3 = {p.x + 0.5 * dt * k2[i].x, p.y + 0.5 * dt * k2[i].y, p.vx + 0.5 * dt * k2[i].vx, p.vy + 0.5 * dt * k2[i].vy};
                
                double x = p_k3.x; double y = p_k3.y; double vx = p_k3.vx; double vy = p_k3.vy;

                double r1 = pow(pow(x+lyap.mu,2)+pow(y,2),1.5); 
                double r2 = pow(pow(x-1+lyap.mu,2)+pow(y,2),1.5); 
                double ax = 2*vy+x-(lyap.mu*(x-1+lyap.mu)/r2)-((1-lyap.mu)*(x+lyap.mu)/r1); 
                double ay = -2*vx+y-(lyap.mu*y/r2)-((1-lyap.mu)*y/r1);

                k3.push_back({p.vx + 0.5 * dt * ax, p.vy + 0.5 * dt * ay, ax, ay});
            }

            for (int i = 0; i < particles.size(); ++i){ // Calculate k4
                Particle& p = particles[i];
                Particle p_k4 = {p.x + dt * k3[i].x, p.y + dt * k3[i].y, p.vx + dt * k3[i].vx, p.vy + dt * k3[i].vy};
                
                double x = p_k4.x; double y = p_k4.y; double vx = p_k4.vx; double vy = p_k4.vy;

                double r1 = pow(pow(x+lyap.mu,2)+pow(y,2),1.5); 
                double r2 = pow(pow(x-1+lyap.mu,2)+pow(y,2),1.5); 
                double ax = 2*vy+x-(lyap.mu*(x-1+lyap.mu)/r2)-((1-lyap.mu)*(x+lyap.mu)/r1); 
                double ay = -2*vx+y-(lyap.mu*y/r2)-((1-lyap.mu)*y/r1);

                k4.push_back({p.vx + dt * ax, p.vy + dt * ay, ax, ay});
            }

            for (int i = 0; i < particles.size(); ++i) {
                Particle& p = particles[i];
                p.x += dt / 6.0 * (k1[i].x + 2.0 * k2[i].x + 2.0 * k3[i].x + k4[i].x);
                p.y += dt / 6.0 * (k1[i].y + 2.0 * k2[i].y + 2.0 * k3[i].y + k4[i].y);
                p.vx += dt / 6.0 * (k1[i].vx + 2.0 * k2[i].vx + 2.0 * k3[i].vx + k4[i].vx);
                p.vy += dt / 6.0 * (k1[i].vy + 2.0 * k2[i].vy + 2.0 * k3[i].vy + k4[i].vy);
            }
        }
};

int main(){

    //===================================== Begin User Input ====================================================
    std::cout << "Reading in user inputs..." << std::endl; std::cout << std::endl; 

    const int NM = 1;                                   // Number of measurements 
    std::string FILE_PATH   = "./Data/Sun-Earth";       // Measurement file path
    std::string FILE_PATH_M = "./Movie Data/Sun-Earth"; // Movie file path

    bool RECORD = true;     // Write MC to .txt file
    bool RECORD_F = true;  // Write frequent MC to .txt file (for trajectories)
    bool MEASURE = true;    // Take discrete measurement updates
    int output_freq = 100;  // Number of steps per output to terminal
    int record_f_freq = 50; // Number of steps per frequent record
    int del_step = 10;      // Number of steps per deletion procedure
    int num_dist = 6;       // Number of distributions recorded per measurement
    //===========================================================================================================

    //===================================== Read in measurement/trajectory info =================================
    std::cout << "Reading in discrete measurements..." << std::endl; std::cout << std::endl; 

    std::ifstream measurement_file(FILE_PATH + "/measurements.txt"); // Measurement file
    double means[NM][DIM]; 
    double stds[NM][DIM];
    int r = 0; 
    std::string line; 
    std::istringstream iss; 
    std::getline(measurement_file, line); // Skip label line
    while((std::getline(measurement_file, line))&&(r < NM)){
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')){
            if(c < DIM){
                means[r][c] = std::stod(line); 
            }else{
                stds[r][c-DIM] = std::stod(line); 
            }
            c++;
        }
        r++;
    }

    Traj lyap; // Trajectory object
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.mu = std::stod(line); // Read in mu
    std::getline(measurement_file, line); // Skip label line
    std::getline(measurement_file, line); lyap.T = std::stod(line);  // Read in T

    std::ifstream size_file(FILE_PATH + "/size.txt"); 
    double size[NM]; r = 0; 
    while(std::getline(size_file, line)) {
        size[r] = std::stod(line);
        r++;
    } 

    MC M; 
    M.epoch = {means[0][0], means[0][1], means[0][2], means[0][3]}; // Grid initial epoch
    M.std = {stds[0][0], stds[0][1], stds[0][2], stds[0][3]};       // Grid width
    M.size = size[0]; 

    double measure_time = lyap.T/NM;            // Time between measurements
    double record_time = measure_time/num_dist; // Time between recording PDF

    measurement_file.close(); size_file.close(); 
    //===========================================================================================================
    std::cout << "Initializing Distribution..." << std::endl; std::cout << std::endl; 

    M.Initialize_MC(lyap);

    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;

    double tt = 0; int nm = 0; int step_count, record_count, record_f_count; double rt, mt; std::ofstream time_file; 
    while(tt < lyap.T){ // Time of uncertainty propagation
        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << tt;
        std::cout << " TU, Active/Total Cells: " << M.size << "/" << M.size << std::endl;
        std::string filename = FILE_PATH + "/M" + std::to_string(nm) + "/MC_0.txt"; M.Record_Data(filename, 0); 
        filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/MC_0.txt"; M.Record_Data(filename, 0); 

        mt = 0; step_count = 0; record_count = 1; record_f_count = 1; 
        while(mt < measure_time){ // Time between measurements 
            rt = 0; std::ifstream time_file(FILE_PATH + "/Times/time" + std::to_string(nm) + ".txt"); 
            double init_time = 0; 

            while(rt < record_time){ // Time between recording the PDF
                std::string t; getline(time_file, t); double dt = std::stod(t)-init_time; init_time = std::stod(t); rt += dt; 

                M.update_prob(lyap, dt); 
                
                step_count+=1; 
                if (step_count % output_freq == 0){ // Print size to terminal
                    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
                    std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt;
                    std::cout << " TU, Active/Total Cells: " << M.size << "/" << M.size << std::endl;
                }    

                if(RECORD_F){ // Record MC Trajectory 
                    if (step_count % record_f_freq == 0){ // Print size to terminal
                        std::string filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/MC_" + std::to_string(record_f_count) + ".txt"; M.Record_Data(filename, mt + rt);
                        record_f_count+=1;
                    }
                }
            }
            if(step_count % output_freq != 0){
                finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
                std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt;
                std::cout << " TU, Active/Total Cells: " << M.size << "/" << M.size << std::endl;
            }

            if(RECORD){ // Record MC Scatter Plot
                std::cout << std::endl; 
                std::cout << "RECORDING PDF AT: " << tt + mt + rt << " TU..." << std::endl;
                std::cout << std::endl; 
                std::string filename = FILE_PATH + "/M" + std::to_string(nm) + "/MC_" + std::to_string(record_count) + ".txt"; M.Record_Data(filename, mt + rt);
                filename = FILE_PATH_M + "/M" + std::to_string(nm) + "/MC_" + std::to_string(record_f_count) + ".txt"; M.Record_Data(filename, mt + rt);
                record_f_count+=1;
                record_count+=1;
            }

            mt += rt; 
        }

        tt += mt; 

        if((MEASURE)&&(tt < lyap.T)){ // Perform discrete measurement update
            std::cout << std::endl; 
            std::cout << "PERFORMING RESAMPLING AT: " << tt << " TU..." << std::endl;
            std::cout << std::endl; 
            nm++; 

            M.size = size[nm];
            M.epoch = {means[nm][0], means[nm][1], means[nm][2], means[nm][3]};
            M.std = {stds[nm][0], stds[nm][1], stds[nm][2], stds[nm][3]}; 
            M.Initialize_MC(lyap);
        }

        time_file.close(); measurement_file.close(); 
    }

    return 0;
}