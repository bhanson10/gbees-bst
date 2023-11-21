/*==============================================================================

MAIN.CPP

==============================================================================*/

#include "GridBayes.h"
#include "GridBayesDerived.h"

#define RECORD_DATA 1 // Record simulation data to text file.

// 1 step = 1 minutes for planetary. Typical LEO orbital period = 1.5 to 2 hrs
int step_record = 400;         
int step_measure = 1200000;   // Measurement every orbital period (approx)
int step_display = 2;

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
  //cout << "Begin list initialization. ";
  auto start = chrono::high_resolution_clock::now();
  
  Lor.initialize();
  
  auto finish = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = finish - start;
  //cout << "Initialization time: " << elapsed.count() << " seconds.\n";
  //cout << "Number of entries in the list: " << Lor.nn - 1 << endl << endl;

  // Record the initial data.
  int count = 0;

  string data_file_name = "./Data/M0/pdf_0.txt";
  Lor.record_data(data_file_name);

  //Lor.measurement_update(); // Initial measurement

  // Start the simulation
  
  //cout << "The simulation begins!" << endl;
  start = chrono::high_resolution_clock::now(); // Start time
  double timer1 = 0;
  double timer2 = 0;

  cout << "Timestep: 0-0, Program time: 0 s, Sim. time: 0";
  cout << " TU, Active/Total Cells: " << Lor.mm << "/" << Lor.nn << endl; 
  for (int i = 1; i <= 2000; i++) {
    Lor.modify_list();      
    Lor.march_PDF();              // March the PDF in time
    Lor.march_RK4(Lor.xs);        // March the actual system in time

    // Measurement update every x steps
    if (i % (step_measure) == 0) {
      Lor.measurement_update();
    }
    
    // Display data to the console every x steps
    if (i % (step_display) == 0) {
      // Printout progress to console
      finish = chrono::high_resolution_clock::now();
      elapsed = finish - start;
      cout << "Timestep: 0-" << i << ", Program time: " << elapsed.count() << " s, Sim. time: " << Lor.t_sim;
      cout << " TU, Active/Total Cells: " << Lor.mm << "/" << Lor.nn << endl; 
      //cout << elapsed.count() << " " << Lor.t_sim << " " << Lor.nn << endl;
    }

    // Record data every x steps
    if (i % (step_record) == 0) {
      // Record every x steps (400/720/43200)
      data_file_name = "./Data/M0/pdf_" + to_string(count+1) + ".txt";
      if (RECORD_DATA) {
        // Record the list data
        Lor.record_data(data_file_name);
      }
      count++;
    }
  }

	return 0;
}

