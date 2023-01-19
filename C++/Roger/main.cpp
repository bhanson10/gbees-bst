/*==============================================================================

MAIN.CPP

==============================================================================*/

#include "GridBayes.h"
#include "GridBayesDerived.h"

#define RECORD_DATA 1 // Record simulation data to text file.

// 1 step = 1 minutes for planetary. Typical LEO orbital period = 1.5 to 2 hrs
int step_record = 400;         
int step_measure = 1200000;   // Measurement every orbital period (approx)
int step_display = 400;

string name_ver = "v7";

double xc_test[6] = { 0,0,0,0,0,0 };
double xlp[6] = { 0,0,0,0,0,0 };
double xlc[6] = { 0,0,0,0,0,0 };

int main() {
  
  // Initialize the class
  // The simulation settings are defined in "GridBayesDerived.h"

  Lorenz3D Lor;
  ofstream simfile;

  Lor.check_init(); // Check if the class has been initialized properly.

  // Array, grid list initialization
  cout << "Begin list initialization. ";
  auto start = chrono::high_resolution_clock::now();
  
  Lor.initialize();
  
  auto finish = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = finish - start;
  cout << "Initialization time: " << elapsed.count() << " seconds.\n";
  cout << "Number of entries in the list: " << Lor.nn - 1 << endl << endl;

  // Record the initial data.
  int count = 0;

  string data_file_name = "cpp_sim_data_" + name_ver + "_" + to_string(count) + ".txt";
  
  Lor.measurement_update(); // Initial measurement

  // Start the simulation
  
  cout << "The simulation begins!" << endl;
  start = chrono::high_resolution_clock::now(); // Start time
  double timer1 = 0;
  double timer2 = 0;

  Lor.modify_list();            // Rearrange list and add points

   if (RECORD_DATA) {
      Lor.record_data(data_file_name);  // Record the initial dataset

      // Simulation states
      simfile.open("simulation_states_" + name_ver + ".txt");
      simfile << Lor.t_sim;
      for (int i = 0; i < Lor.N_DIM; i++) {
          simfile << "," << Lor.xs[i];
      }
      simfile << endl;
  }

  for (int i = 1; i <= 2000; i++) {
    Lor.march_PDF();              // March the PDF in time
    Lor.march_RK4(Lor.xs);        // March the actual system in time

    // Measurement update every x steps
    if (i % (step_measure) == 0) {
      Lor.measurement_update();
    }
    Lor.modify_list();            // Rearrange list and add points
    
    

    // Display data to the console every x steps
    if (i % (step_display) == 0) {
      // Printout progress to console
      finish = chrono::high_resolution_clock::now();
      elapsed = finish - start;
      cout << "Step: " << i;
      cout << ", program time: " << elapsed.count() << " seconds";
      cout << ", sim time: " << Lor.t_sim << " seconds";
      cout << ", active/total cells: " << Lor.mm << " / " << Lor.nn << ".";
      cout << endl;

    }

    // Record data every x steps
    if (i % (step_record) == 0) {
      // Record every x steps (400/720/43200)
      count++;
      data_file_name = "cpp_sim_data_" + name_ver + "_" + to_string(count) + ".txt";
      if (RECORD_DATA) {
        // Record the list data
        Lor.record_data(data_file_name);

        // Record simulation planetary states
        simfile << Lor.t_sim;
        for (int i = 0; i < Lor.N_DIM; i++) {
          simfile << "," << Lor.xs[i];
        }
        simfile << endl;
      }
    }
  }

  // Record every x steps (400/720/43200)
  count++;
  data_file_name = "cpp_sim_data_" + name_ver + "_" + to_string(count) + ".txt";
  if (RECORD_DATA) {
    // Record the list data
    Lor.record_data(data_file_name);

    // Record simulation planetary states
    simfile << Lor.t_sim;
    for (int i = 0; i < Lor.N_DIM; i++) {
      simfile << "," << Lor.xs[i];
    }
    simfile << endl;
  }

  if (RECORD_DATA) simfile.close();

	return 0;
}

