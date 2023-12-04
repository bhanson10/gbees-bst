# Grid-based, Bayesian Estimation Exploiting Sparsity (GBEES)
The "GBEES" repsository containes the codebase accompying the paper "State Estimation of Chaotic Trajectories: A Higher-Dimensional, Grid-Based, Bayesian Approach to Uncertainty Propagation" presented at the January 2024 AIAA/AAS Space Flight Mechanics Meeting. Below is an in-depth summary explaining the proper usage of all the software included in said repository. In all, the codebase provides the computational framework necessary to: <br> 
1. Efficiently, continuous-time propagate an initially Gaussian distribution governed by the dynamics of a nonlinear system <br>
2. Perform discrete measurement updates at given epochs via Bayes theorem <br>
3. Compare the accuracy of said propagation with a Monte Carlo (MC) simulation of similar resolution <br>
4. Compare the efficiency of the GBEES propagation with the MC simulation of similar resolution  <br> 

The current framework is applicable to 4D systems, and is specifically utilized to propagate orbital uncertainty in the PCR3BP. Users can use any initial conditions, period, and celestial system as long as the data files are placed in the correct locations (more on this later). 

## General Formulation
There are two provided systems that GBEES is applied to: <br> 
1. 3D Lorenz Chaotic Attractor <br>
2. 4D Planar Circular Restricted Three Body Problem <br>

Each of these systems has their own folder within the GBEES repository. GBEES is propagated via the C++ code _main.cpp_ in the "/GBEES" subfolder within each of these folder. An MC simulation of similar resolution as the GBEES propagation is propagated via the C++ code _main.cpp_ in the "/MC" subfolder. These results are then compared using the _plot_BLANK_m._ MATLAB code. The MATLAB code creates all of the figures from the paper. The _movie_BLANK_m._ stitches together a series of images to create an animation of the PDFs spreading over phase space, and can be created by saving the PDF data at a higher frequency than is defaultly done. <br> 

Here is a breakdown of all of the subrepositories within each of the three examples provided, as well as their functions: <br>

* '/GBEES/Data': This is where the PDFs that are written to _.txt_ files during GBEES are stored. Each PDF is sequestered into its specific measurement number subrepository, '/M#'. Ensure that the expected number of measurements taken are equal to the number of '/M#' subfolders. The discrete measurement updates are also saved in this file, in the form of a _measurement_BLANK.txt_ file. Update this .txt file with the measurement vector. Additionally, the size of the timesteps taken by the explicit scheme are saved in this folder, under the '/Data/Times' repository. This information is used by the MC simulation, so ensure once your GBEES simulation has run, you copy the times over to the '/MC/Times' file. <br>
* '/GBEES/Movie Data': This is folder has a similar structure to the '/GBEES/Data', but with a more frequent saving procedure of the PDFs. This way, the PDFs can be stitched together via MATLAB to create an animation of the PDF expanding and moving over time. To adjust the frequency of saving the PDF, change the _record_time_ variable in the '/GBEES/main.cpp' file so that it is lower, but still a factor of _measure_time_. Additionally, ensure that the _Record_Data_ function is writing to the correct location by changing the folder to './Movie Data'. <br>
* '/MC/Data': This folder has a similar function as the '/GBEES/Data' folder, but for the MC simulation. Rather than storing a PDF, all of the particles and their current locations at a given epoch are stored in '/Data/M#', again where 'M#' is based on what measurement number the simulation is on. The black dots in the figures from the paper that compare the PDFs and the MC simulation come from this folder. <br>
* '/MC/Movie Data': This folder has a similar function as the '/MC/Data' folder, but stores the particle locations in phase space at a more frequent rate. This way, the full trajectory may be stored and added to the MC comparison figures, which are the gray lines in the paper. <br>
* '/Initial Conditions': This folder in thePCR3BP examples contain _.csv_ files pulled from the JPL Three-Body Periodic Orbit catalog which have the initial conditions of the Lyapunov trajectories. <br>

Here is a breakdown of all of the files within the examples provided, as well as their usages: <br>
* _measurements.txt_: This _.txt_ file contains the discrete measurement vectors that are fed to the GBEES and MC simulations at the provided epoch. They are of the form 'x y z' with spaces in between float values. Each line represents a new measurement.
* _size.txt_: This _.txt_ file contains the number of particles that should be initialized by the MC simulation at each discrete measurement update and is equal to the number of active cells in the final distributions of the GBEES simulation, prior to discrete measurement update. This way, the GBEES and MC simulations have equal resolution, thus can be compared for efficiency. 
* _runtime.txt_: This _.txt_ file stores the information necessary to compare the runtimes and sizes of the GBEES and MC simulations (more on this later). 
Now that each of the folders and files within the GBEES repositiory has been explained, we'll go into the intracies of each provided example, as they all have slight differences that are implemented. <br>

## '/Lorenz3D'
The 3D Lorenz attractor example provides two different implementations: <br> 
1. The legacy implementation, stored in the '/Legacy' subfolder, which was the implementation utilized in the GBEES paper by [T. Bewley of 2012](https://www.sciencedirect.com/science/article/pii/S0005109812000908). <br>
2. The current implementation, stored in the '/Current' subfolder, which includes all the modifications made to improve the efficiency, discussed in the associated AIAA/AAS conference paper. <br>

These two implementations are compared via the '/Timing' folder (more on this later). 

## '/Jupiter-Europa'
This example takes a set of initial conditions that result in a Lyapunov orbit about the L3 Jupiter-Europa libration point, sourced from the JPL Three-Body Periodic Orbit catalog. We propagate using GBEES for an entire period, corresponding to about 3.5 days, then compare this propagation with a MC simulation. For this example, we assume no epistemic uncertainty, thus the diffusion term is 0, and only advection is considered. 

## '/Sun-Earth'
This example takes a set of initial conditions that result in a Lyapunov orbit about the L3 Sun-Earth libration point, sourced from the JPL Three-Body Periodic Orbit catalog. For this simulation, we propagate using GBEES through the close approach, corresponding to about 1.8 days, then compare this propagation with a MC simulation. For this example, we assume epistemic uncertainty, thus including the diffusion term. 

## '/Timing'
This folder is dedicated to comparing the results of GBEES and the corresponding MC simulation. To do this, copy and paste the output of _main.cpp_ from both the GBEES and MC simulations and store them in the _runtime.txt_ files. Then, copy these files over to the '/Timing' folder, and save them in the corresponding subfolder (or make your own, for a new simulation). Using the _compare_times.m_ MATLAB file, edit the path of the folder that you would like to compare, as well as the axes units, and then run the program. This will output the efficiency comparison figures from the paper. 
<br><br>
For further information about code usage, please contact blhanson@ucsd.edu
