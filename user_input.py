## File that user has to provide for FSSH simulations performed in FSSH.py
## Contains parameters for simulation, potential and the initial conditions

import numpy as np

nquant=2 ## Quantum states
nclass=1 ## number of classical d.o.f.
ntraj=100 ## number of FSSH trajectories

## FSSH state of the system variables
mass=np.zeros(nclass)
x=np.zeros(nclass)
v=np.zeros(nclass)
ci=np.zeros(nquant,dtype=complex)

## Parameters of potential
mass[0]=2000
hbar=1.0
AA=0.01
BB=0.6
CC=0.001
DD=1.0
E0=0.05

## Parameters of simulation
dtc=20.0        ## Classical time step
dtq=1.0         ## Quantum time step
total_time=2000.0   ## Total simulation tiem
ntim=int(total_time/dtc)    ## Number of steps

## Set nprint to 1 to print detailed information for each trajectory
nprint=0

###############################################
## Input: x(nclass)
## Returns: Potential and its derivative
## Potential V[nquant,nquant], derivative: dv_dx[nquant,nquant,nclass]
def pot(x):
    V=np.zeros((nquant,nquant))
    dv_dx=np.zeros((nquant,nquant,nclass))
    V[0,0]=AA*np.tanh(BB*x)
    V[1,1]=-V[0,0]
    #V[1,1]=AA/2*np.tanh(BB*x)
    V[0,1]=CC*np.exp(-DD*x*x)
    V[1,0]=V[0,1]

    dv_dx[0,0,0]=AA*BB*(1-np.tanh(BB*x)**2)
    dv_dx[1,1,0]=-dv_dx[0,0,0]
    dv_dx[0,1,0]=-2*CC*DD*x*np.exp(-DD*x*x)
    dv_dx[1,0,0]=-dv_dx[0,1,0]

    return(V,dv_dx)
###################################################

## Initial condition
## Returns position, velocity, initial state, and ci's
## Choosing position from a Gaussian distribution exp(-(x-x0)/(2sigma_x**2)), 
## x0=-5 au, sigma_x=1/sqrt(2) au, consistent with the corresponding exact QM calculation
## momentum is chosen with a distribution with sigma_x.sigma_p=hbar/2
def init_cond():
    k=np.sqrt(2*mass*0.03)
    sig_x=1.0/np.sqrt(2.0)
    sig_v=0.5/mass[0]*hbar/sig_x
    x[0]=-5.0 + sig_x*np.random.normal()
    v[0]=hbar*k/mass[0] + sig_v*np.random.normal()
    lamda=0
    ci[:]=0
    ci[lamda]=1.0
    return(x,v,lamda,ci)
###################################################
