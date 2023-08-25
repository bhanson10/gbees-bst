/*==============================================================================

MAIN.CPP

==============================================================================*/

#include "GridBayes.h"
#include "GridBayesDerived.h"

#define RECORD_DATA 1 // Record simulation data to text file.

// 1 step = 1 minutes for planetary. Typical LEO orbital period = 1.5 to 2 hrs
int step_record = 1;         
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
  
  Lor.measurement_update(); // Initial measurement

  // Start the simulation
  
  cout << "The simulation begins!" << endl;
  start = chrono::high_resolution_clock::now(); // Start time

  Lor.modify_list();            // Rearrange list and add points

  std::string filename = "time_complexity_nl.txt"; std::ofstream myfile; myfile.open(filename);
  for (int i = 1; i <= 2000; i++) {
    Lor.march_PDF();              // March the PDF in time
    Lor.march_RK4(Lor.xs);        // March the actual system in time

    Lor.modify_list();            // Rearrange list and add points

    // Record data every x steps
    if (i % (step_record) == 0) {
        // Record every x steps (400/720/43200)
        auto finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
        myfile << Lor.nn << " " << elapsed.count() << std::endl;
        start = finish; 
    }
  }

	return 0;
}

