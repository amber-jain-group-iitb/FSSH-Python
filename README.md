Python codes to run a basic FSSH simulation. Contains FSSH.py and user_input.py. Original reference: Tully, JCP 93, 1061, 1990

Follows the algorithm provided in Section 9 in the ACS omega review by Jain and Sindhu.
Note that these codes do NOT include decoherence, correct treatment of frustrated hops, or correct calculation of electronic density matrix. Codes with decoherence (from augmented FSSH approach) and correct treatment of frustrated hops and density matrix can be found at https://github.com/amber-jain-group-iitb/AFSSH

Results of this simulations can be compared to numerically exact results from https://github.com/amber-jain-group-iitb/Exact-QD-Python

Running the codes:
  Download both FSSH.py and user_input.py in the same location.
  Requires Python3, with numpy, scipy, and matplotlib libraries. All are freely available.
  Running on terminal: cd to the directory where the above 2 files are stored. Issue on the terminal: python3 FSSH.py
  Can use Python editors (such as Pycharm) to run the code if not comfortable with terminal.
  The codes will reproduce the results for the potential in Fig. 1(a) (non-parallel surfaces) of the review. The code will run 100 trajectories with dtc=20 a.u., and total time=2000 a.u.
  On running, the code will print the initial position and velocity for each trajectory, and after 100 trajectories will create a plot for the absolute values of the density matrix as a function of time (Fig. 3 of review).
  The code should take 1-2 mins on a regular desktop.

Structure of the code:
  FSSH.py - contains the main code for the FSSH algorithm. Detailed algorithm in the review (Sec. 9).
  user_input.py - contains all required input from the user. This includes nquant (number of quantum state), nclass (number of classical d.o.f.), ntraj (number of trajs.), mass, time steps (dtc - classical, dtq - quantum), total simulation time.
    pot(x) in user_input.py takes the classical x (position) as input and returns the potential V[nquant,nquant] and its derivative dv_x[nquant,nquant,nclass]
    init_cond() - this is called for each trajectory. returns the initial FSSH state of the system: x, v (position, velocity) ci (quantum coefficients) and lamda (active state).
  Setting variable nprint=1 will print detailed information for each trajectory. Useful for debugging (for ex. can see if energy is conserved or not).   
  Typically values: ntraj~100, dtc~20 au, dtq~1 au - though these values can depend significantly on the system under study. Make sure energy is conserved and final results are converged w.r.t. the parameters.
 
Things to try (excercies):
  Change dtc and dtq and re-run to check if the plot of density matrix changes.
  Current code is provided for the potential V[0,0]=A tanh(Bx), V[1,1]=-V[0,0], V[0,1]=V[1,0]=C.exp(-D.x^2) (x is the classical position).
  This can be found in the pot(x) function in user_input.py
  Change the potential to V[1,1]=V[0,0]/2. CHANGE the corresponding derivative dv_dx[1,1,0] as well. Re-run the code to reproduce Fig. 4.
  Change the initial energy. Currently the code starts with kinetic energy of 0.03 a.u. Find the variable k in init_cond() in user_input.py and change initial KE to 0.02 a.u and see the differences
  
  Change potential to the 3 Tully models given in Tully, JCP 93, 1061, 1990 to reproduce the original results.
  
  Adding decoherence to the code: (will require basic Python coding)
    Go to function decoherence() in FSSH.py. This is being called at every classical time step. Currently this is not doing anything. Using your favorite decoherence scheme, calculate the probability of a collapse: p_coll. Call a uniform random number rnd between 0 and 1. If rnd< p_coll, set ci[lamda]=1 and ci[i]=0 for i not equal to lamda.
  
  Reversing velocity on frustrated hop: (will require basic Python coding)
    Go to perform_hop() in FSSH.py. Find the condition for frustrated hop within this function. Setting gama=bb/aa (instead of 0) will always reverse velocity along the dij direction.
    On a frustrated hop, calculate the forces of the active surface and the surface to which hop was attempted to. If these forces have opposite signs, set gama=bb/aa (reverse velocity), else gama=0 (do not reverse velocity).
    On reversing velocity on frustrated hop: Jasper and Truhlar, CPL 369, 60, 2003.
  
  Calculating electronic density matrix: (will require basic Python coding)
    Go to compute_averages() in FSSH.py. Here rho[i,j] is calculated as conjg(ci[i])*ci[j].
    An alternate way to calculate rho is: rho[lamda,lamda]=1, rho[i,i]=0 if i.ne.lamda, and rho[i,j]=conjg(ci[i])*ci[j] if i.ne.j
    This way to compute populations from Landry, Falk, Subotnik, JCP 139, 211101, 2013.
