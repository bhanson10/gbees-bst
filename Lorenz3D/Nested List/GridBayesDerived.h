#ifndef GRIDBAYESDERIVED__H
#define GRIDBAYESDERIVED__H

#include "GridBayes.h"

/*==============================================================================
  DERIVED CLASS DEFINITIONS
==============================================================================*/

// 3 states Lorenz equation ----------------------------------------------------

class Lorenz3D : public GridBayes {

public:
  // Three states Lorenz equation

  // Constructors	and destructors, initialize const values
  Lorenz3D(
    // Base class constant values.
    int n_pmax = 100000
    , int n_upd = 1
    , double t_max = 1
    , double dt = 5e-4
    , double pdf_thresh = 2e-5
    , vector<double> dx = vector<double>(3, 0.5)
    , vector<double> x0 = { -11.5, -10, 9.5 }
    , vector<double> x0_stdev = { 0.5, 0.5, 0.5 }
    , int dn = 5
  )
    : GridBayes(3, n_pmax, n_upd, t_max, dt, pdf_thresh, dx, x0, x0_stdev, dn)
  {}

  // Constant variables (Lorenz model parameters)
  double SIGMA = 4;
  double B = 1;
  double R = 48;

  // Virtual functions
  void model_f(double* f, const double* x, const double* xh)
  {
    // dx/dt = f(x,t) model : 3 states Lorenz equation
    f[0] = SIGMA * (x[1] - (x[0] + xh[0]));
    f[1] = -(x[1] + xh[1]) - x[0] * x[2];
    f[2] = -B * (x[2] + xh[2]) + x[0] * x[1] - B * R;
  }
  double measurement_model_pdf(int id)
  {
    // Measurement model 
    // Returns the measurement pdf for the list member [id]
    
    // Measure z axis only
    double z_meas = xs[2]; // simulated z state value
    double stdev = 0.5;   // measurement standard deviation

    // Current position in the list member [id]
    double z_list = list.pos[id][2] * DX[2];

    // Determine the pdf value based on Gaussian white noise model.
    double pdf = exp( -pow(z_list - z_meas, 2) / 2.0 / stdev / stdev );

    return pdf;
  }
  double flux_limiter(double x)
  {
    // Monotonized Central Difference (MC)
    // min([ (1+x)/2, 2, 2*x])
    double min = (1 + x) / 2.0;
    if (min > 2) min = 2;
    if (min > 2 * x) min = 2 * x;
    if (min > 0) return min;
    else return 0;

    // van Leer
    // return (x + abs(x)) / (1.0 + abs(x)); 
  }

private:

};

#endif // __GRIDBAYESDERIVED__H