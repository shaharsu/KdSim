This Matlab script package contains all the tools needed to simulate protein-protein interactions
detected by FRET and induced by cellular volume changes, as specified in the original paper.
The package contains the following scripts:

KdSimGUI.m - GUI to help visualize simulation results
run_simpleKd.m - Main function used to generate simulations
simpleKd_step.m - integrator used in conjunction with run_simpleKd.m
minimizeRun_simpleKd.m - Function to iteratively run run_simpleKd.m, varying input to fit experimental results
sseRun_simpleKd.m - minimization function to be used with minimizeRun_simpleKd.m
KdSimGUI.fig - GUIDE file for GUI

If you have any questions or problems with this program, please contact Shahar Sukenik, shaharsu@gmail.com