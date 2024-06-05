# Grid-based, Bayesian Estimation Exploiting Sparsity (GBEES)
The "GBEES" repository contains the codebase accompanying the paper "State Estimation of Chaotic Trajectories: A Higher-Dimensional, Grid-Based, Bayesian Approach to Uncertainty Propagation" presented at the January 2024 AIAA/AAS Space Flight Mechanics Meeting. Below is an in-depth summary explaining the proper usage of all the software included in said repository. In all, the codebase provides the computational framework necessary to: <br> 
1. Efficiently, continuous-time propagate an initially Gaussian distribution governed by the dynamics of a nonlinear system <br>
2. Perform discrete measurement updates at given epochs via Bayes theorem <br>
3. Compare the accuracy of said propagation with a Particle Filter (PF)/Monte Carlo (MC) simulation of similar resolution <br>
4. Compare the efficiency of the GBEES propagation with the PF/MC simulation of high-resolution  <br> 

The current framework is applicable to 4D systems, and is specifically utilized to propagate orbital uncertainty in the PCR3BP. Users can use any initial conditions, period, and celestial system as long as the data files are placed in the correct locations (more on this later). 

## General Formulation
There are three provided systems that GBEES is applied to: <br> 
1. 3D Lorenz Chaotic Attractor <br>
2. 4D Planar Circular Restricted Three Body Problem <br>
3. 6D Restricted Two-Body (Keplerian) Problem <br>

Each of these systems has their own folder within the GBEES repository. GBEES is propagated via the C++ code _main.cpp_ in the "/GBEES" subfolder within each of these folder. A PF/MC simulation of similar resolution as the GBEES propagation is propagated via the C++ code _main.cpp_ in the "/PF" or "/MC" subfolder. <br>

Here is a breakdown of all of the subrepositories within each of the three examples provided, as well as their functions: <br>

* '/GBEES/Epochs': This is where the PDFs that are written to _.txt_ files during GBEES are stored. Each PDF is sequestered into its specific simulation number subrepository, '/M#'.The discrete measurement updates are also saved in this file, in the form of a _measurement.txt_ file. Update this .txt file with the measurement vector. Additionally, the size and the number of timesteps taken by the explicit scheme are saved in this folder, in the _timing.txt_ file. This information is used for comparing computational efficiency. <br>


Here is a breakdown of all of the files within the examples provided, as well as their usages: <br>
* _measurements.txt_: This _.txt_ file contains the discrete measurement vectors that are fed to the GBEES and MC simulations at the provided epoch, the standard deviation associated with each measurement, the gravitational coefficient of the system, and the period of the orbit (or the time of total propagation). This should be the form of the measurement file (for the PCR3BP): <br> <br>

<div align="center">
  
|           x (LU)           |           y (LU)           |          vx (LU/TU)           |          vy (LU/TU)           | 
|:--------------------------:|:--------------------------:|:-----------------------------:|:-----------------------------:|
|             x0             |             y0             |              vx0              |              vy0              |
|                            |                            |                               |                               |
|          dx (LU)           |          dy (LU)           |          dvx (LU/TU)          |          dvy (LU/TU)          | 
|        $\sigma_x^2$        |  $\sigma_x\cdot\sigma_y$   |  $\sigma_x\cdot\sigma_{vx}$   |  $\sigma_x\cdot\sigma_{vy}$   |
|  $\sigma_y\cdot\sigma_x$   |        $\sigma_y^2$        |  $\sigma_y\cdot\sigma_{vx}$   |  $\sigma_y\cdot\sigma_{vy}$   |
| $\sigma_{vx}\cdot\sigma_x$ | $\sigma_{vx}\cdot\sigma_y$ |        $\sigma_{vx}^2$        | $\sigma_{vx}\cdot\sigma_{vy}$ |
| $\sigma_{vy}\cdot\sigma_x$ | $\sigma_{vy}\cdot\sigma_y$ | $\sigma_{vy}\cdot\sigma_{vx}$ |        $\sigma_{vy}^2$        |
|                            |                            |                               |                               |
|           $\mu$            |                            |                               |                               |
|          $\mu_1$           |                            |                               |                               |
|           T (TU)           |                            |                               |                               |
|             T              |                            |                               |                               |

</div>

* _timing.txt_: This _.txt_ file stores the information necessary to compare the runtimes and sizes of the GBEES and PF/MC simulations (more on this later). 
Now that each of the folders and files within the GBEES repository has been explained, we'll go into the intricacies of each provided example, as they all have slight differences that are implemented. <br>

## '/Lorenz3D'
The example uses the 3D Lorenz attractor dynamics, and includes all the modifications made to improve the efficiency, discussed in the associated AIAA/AAS conference paper. <br>


## '/R2BP'
This example uses the Restricted Two-Body (Keplerian) dynamics, using an initial condition that results in an elliptical orbit about a main body (either Earth or Europa) and propagate for multiple revolutions without measurement updates until the uncertainty becomes highly non-Gaussian. For this simulation, we compare a UKF (Cartesian and equinoctial) and EnKF (Cartesian and equinoctial) to a high-resolution PF (truth model).

## '/PCR3BP'
This example uses the Planar Circular Restricted Three-Body dynamics, using an initial condition that results in a Low-Prograde Orbit about Europa, sourced from the JPL Three-Body Periodic Orbit catalog. For this simulation, we propagate using GBEES until positional uncertainty becomes non-Gaussian, corresponding to about 14 hours, then compare this propagation with a UKF, EnKF, and a high-resolution PF (truth model). 

<br><br>
For further information about code usage, please contact blhanson@ucsd.edu
