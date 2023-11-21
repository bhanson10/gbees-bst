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
const int DIM = 3; 
const int NM = 1; 

struct Particle {
    double x, y, z;
};

class Traj{
    public: 
        double sigma; 
        double b;
        double r;
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

        void Initialize_MC(Traj Lor){
            particles.clear(); 

            std::default_random_engine generator;
            std::normal_distribution<double> distribution1(epoch[0],std[0]);
            std::normal_distribution<double> distribution2(epoch[1],std[1]);
            std::normal_distribution<double> distribution3(epoch[2],std[2]);

            for(int q = 0; q < size; q++){
                Particle p;
                p.x = distribution1(generator);
                p.y = distribution2(generator);
                p.z = distribution3(generator);

                particles.push_back(p); 
            }
        }

        void Record_Data(std::string filename, double t){
            std::ofstream myfile; myfile.open(filename);
            myfile << t << std::endl; 

            for (Particle& p : particles){
                myfile << p.x << " " << p.y << " " << p.z << " " << std::endl; 
            }
        }

        void update_prob(Traj Lor, double dt){
            std::vector<Particle> k1, k2, k3, k4;

            for (Particle& p : particles){ // Calculate k1
                double x = p.x; double y = p.y; double z = p.z;

                k1.push_back({Lor.sigma*(y-x), -y-x*z, -Lor.b*z+x*y-Lor.b*Lor.r});
            }

            for (int i = 0; i < particles.size(); ++i){ // Calculate k2
                Particle& p = particles[i];
                Particle p_k2 = {p.x + 0.5 * dt * k1[i].x, p.y + 0.5 * dt * k1[i].y, p.z + 0.5 * dt * k1[i].z};
                
                double x = p_k2.x; double y = p_k2.y; double z = p_k2.z;

                k2.push_back({Lor.sigma*(y-x), -y-x*z, -Lor.b*z+x*y-Lor.b*Lor.r});
            }

            for (int i = 0; i < particles.size(); ++i){ // Calculate k3
                Particle& p = particles[i];
                Particle p_k3 = {p.x + 0.5 * dt * k2[i].x, p.y + 0.5 * dt * k2[i].y, p.z + 0.5 * dt * k2[i].z};

                double x = p_k3.x; double y = p_k3.y; double z = p_k3.z;

                k3.push_back({Lor.sigma*(y-x), -y-x*z, -Lor.b*z+x*y-Lor.b*Lor.r});
            }

            for (int i = 0; i < particles.size(); ++i){ // Calculate k4
                Particle& p = particles[i];
                Particle p_k4 = {p.x + dt * k3[i].x, p.y + dt * k3[i].y, p.z + dt * k3[i].z};

                double x = p_k4.x; double y = p_k4.y; double z = p_k4.z;

                k4.push_back({Lor.sigma*(y-x), -y-x*z, -Lor.b*z+x*y-Lor.b*Lor.r});
            }

            for (int i = 0; i < particles.size(); ++i) {
                Particle& p = particles[i];
                p.x += dt / 6.0 * (k1[i].x + 2.0 * k2[i].x + 2.0 * k3[i].x + k4[i].x);
                p.y += dt / 6.0 * (k1[i].y + 2.0 * k2[i].y + 2.0 * k3[i].y + k4[i].y);
                p.z += dt / 6.0 * (k1[i].z + 2.0 * k2[i].z + 2.0 * k3[i].z + k4[i].z);
            }
        }
};

int main(){
    std::ifstream measurement_file("measurements.txt");
    std::ifstream size_file("size.txt"); 
    double measurements[NM][DIM]; int r = 0; std::string line; std::istringstream iss;
    while(std::getline(measurement_file, line)) {
        int c = 0;
        std::istringstream iss(line);
        while (std::getline(iss, line, ' ')) {
            measurements[r][c] = std::stod(line); 
            c++;
        }
        r++;
    }
    double size[NM]; r = 0; 
    while(std::getline(size_file, line)) {
        size[r] = std::stod(line);
        r++;
    } 
    //===================================== Begin User Input ======================================
    MC M; M.std  = {1, 1, 1}; M.size = size[0]; 
    M.epoch = {measurements[0][0], measurements[0][1], measurements[0][2]};
    Traj Lor; Lor.sigma = 4; Lor.b = 1; Lor.r = 48; Lor.T = 1; 
    double measure_time = Lor.T; double record_time = measure_time/5; bool RECORD = true; bool RECORD_F = true; bool MEASURE = true;
    int output_freq = 100; int record_f_freq = 5;
    //====================================== End User Input =======================================
    auto start = std::chrono::high_resolution_clock::now(); 
    std::chrono::duration<double> elapsed;

    M.Initialize_MC(Lor);

    double tt = 0; int nm = 0; int record_count; int record_f_count; double rt; double mt; int step_count; std::ofstream time_file; 
    while(tt < Lor.T){ // Time of uncertainty propagation
        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
        std::cout << "Timestep: " << nm << "-0, Program time: " << elapsed.count() << " s, Sim. time: " << tt << " TU" << std::endl;
        std::string filename = "./MC Data/M" + std::to_string(nm) + "/MC_0.txt"; M.Record_Data(filename, 0); 
        filename = "./MC Freq Data/M" + std::to_string(nm) + "/MC_0.txt"; M.Record_Data(filename, 0); 

        mt = 0; step_count = 0; record_count = 1; record_f_count = 1; 
        while(mt < measure_time){ // Time between measurements 
            rt = 0; std::ifstream time_file("./Times/time" + std::to_string(nm) + ".txt"); 
            double init_time = 0; 

            while(rt < record_time){ // Time between recording the PDF
                std::string t; getline(time_file, t); double dt = std::stod(t)-init_time; init_time = std::stod(t); rt += dt; 

                M.update_prob(Lor, dt); 
                
                step_count+=1; 
                if (step_count % output_freq == 0){ // Print size to terminal
                    auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
                    std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt << " TU" << std::endl; 
                }    

                if(RECORD_F){ // Record MC Trajectory 
                    if (step_count % record_f_freq == 0){ // Print size to terminal
                        std::string filename = "./MC Freq Data/M" + std::to_string(nm) + "/MC_" + std::to_string(record_f_count) + ".txt"; M.Record_Data(filename, mt + rt);
                        record_f_count+=1;
                    }
                }
            }
            if(step_count % output_freq != 0){
                finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start;
                std::cout << "Timestep: " << nm << "-" << step_count << ", Program time: " << elapsed.count() << " s, Sim. time: " << tt + mt + rt << " TU" << std::endl; 
            }
 
            if(RECORD){ // Record MC 
                //std::cout << std::endl; 
                //std::cout << "RECORDING MC AT: " << tt + mt + rt << " TU..." << std::endl;
                //std::cout << std::endl; 
                std::string filename = "./MC Data/M" + std::to_string(nm) + "/MC_" + std::to_string(record_count) + ".txt"; M.Record_Data(filename, mt + rt);
                filename = "./MC Freq Data/M" + std::to_string(nm) + "/MC_" + std::to_string(record_f_count) + ".txt"; M.Record_Data(filename, mt + rt);
                record_f_count+=1;
                record_count+=1;
            }

            mt += rt; 
        }

        tt += mt; 

        if((MEASURE)&&(tt < Lor.T)){ // Perform discrete measurement update
            std::cout << std::endl; 
            std::cout << "PERFORMING BAYESIAN UPDATE AT: " << tt << " TU..." << std::endl;
            std::cout << std::endl; 
            nm++; 

            M.size = size[nm];
            M.epoch = {measurements[nm][0], measurements[nm][1], measurements[nm][2]};
            M.Initialize_MC(Lor);
        }

        time_file.close(); measurement_file.close(); 
    }

    return 0;
}