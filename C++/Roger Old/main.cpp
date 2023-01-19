/*==============================================================================

MAIN.CPP

==============================================================================*/

#include "GridBayes.h"
#include "GridBayesDerived.h"
#include <fstream>

// 1 step = 1 minutes for planetary. Typical LEO orbital period = 1.5 to 2 hrs
int step_record = 30;         
int step_measure = 1200000;   // Measurement every orbital period (approx)
int step_display = 30;

double xc_test[6] = { 0,0,0,0,0,0 };
double xlp[6] = { 0,0,0,0,0,0 };
double xlc[6] = { 0,0,0,0,0,0 };

int main() {

  // Initialize the class
  // The simulation settings are defined in "GridBayesDerived.h"

  Lorenz3D Lor;
  
  Lor.check_init(); // Check if the class has been initialized properly.

  // Array, grid list initialization
  cout << "Begin list initialization. ";
  auto start = chrono::high_resolution_clock::now();
  
  Lor.initialize();
  cout << endl << "Number of Entries in PDF: " << sizeof(Lor.list.pdf) << endl;
  auto finish = chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed = finish - start;
  cout << "Initialization time: " << elapsed.count() << " seconds.\n";
  cout << "Number of entries in the list: " << Lor.nn - 1 << endl << endl;
  
  Lor.measurement_update(); // Initial measurement

  // Start the simulation
  
  cout << "The simulation begins!" << endl;
  start = chrono::high_resolution_clock::now(); // Start time
  double timer1 = 0;
  double timer2 = 0;

  Lor.modify_list();            // Rearrange list and add points
  

  ofstream pdffile("gbees_cpp_pdf.txt");
  ofstream posfile_i("gbees_cpp_pos_i.txt");
  ofstream posfile_j("gbees_cpp_pos_j.txt");
  ofstream posfile_k("gbees_cpp_pos_k.txt");
  

  for(int i = 0; i < sizeof(Lor.list.pdf); i++){
        pdffile << Lor.list.pdf[i] << " " ;
  }
  pdffile << endl;

  for(int i = 0; i < sizeof(Lor.list.pos); i++){
        posfile_i << Lor.list.pos[i][0] << " ";
  }
  posfile_i << endl;

  for(int i = 0; i < sizeof(Lor.list.pos); i++){
        posfile_j << Lor.list.pos[i][1] << " ";
  }
  posfile_j << endl;

  for(int i = 0; i < sizeof(Lor.list.pos); i++){
        posfile_k << Lor.list.pos[i][2] << " ";
  }
  posfile_k << endl;

  //  for (int i = 1; i <= Lor.T_MAX / Lor.DT; i++) {

  for (int i = 1; i <= Lor.T_MAX / Lor.DT; i++) {
      //Lor.modify_list();            // Rearrange list and add points
      cout << endl << "Number of Entries in PDF: " << sizeof(Lor.list.pdf) << endl;
      Lor.march_PDF();              // March the PDF in time
      Lor.march_RK4(Lor.xs);        // March the actual system in time

      Lor.measurement_update();

      if(i % 400 == 0){
        for(int j = 0; j < sizeof(Lor.list.pdf); j++){
          pdffile << Lor.list.pdf[j] << " " ;
        }
        pdffile << endl;

        for(int j = 0; j < sizeof(Lor.list.pos); j++){
          posfile_i << Lor.list.pos[j][0] << " ";
        }
        posfile_i << endl;

        for(int j = 0; j < sizeof(Lor.list.pos); j++){
          posfile_j << Lor.list.pos[j][1] << " ";
        }
        posfile_j << endl;

        for(int j = 0; j < sizeof(Lor.list.pos); j++){
          posfile_k << Lor.list.pos[j][2] << " ";
        }
        posfile_k << endl;
      }
  }

  Lor.modify_list();            // Rearrange list and add points

  for(int j = 0; j < sizeof(Lor.list.pdf); j++){
    pdffile << Lor.list.pdf[j] << " " ;
  }
  pdffile << endl;

  for(int j = 0; j < sizeof(Lor.list.pos); j++){
    posfile_i << Lor.list.pos[j][0] << " ";
  }
  posfile_i << endl;

  for(int j = 0; j < sizeof(Lor.list.pos); j++){
    posfile_j << Lor.list.pos[j][1] << " ";
  }
  posfile_j << endl;

  for(int j = 0; j < sizeof(Lor.list.pos); j++){
    posfile_k << Lor.list.pos[j][2] << " ";
  }
  posfile_k << endl;

  return 0;
}

