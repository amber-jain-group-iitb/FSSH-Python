## Surface hopping code
print("WARNING: DECOHERENCE IS NOT INCLUDED")
print("WARNING: ON FRUSTRATED HOP, STATE WILL NOT CHANGE AND VELOCITY WILL NOT BE REVERSED")
print("WARNING: Calculating electronic density matrix as rho(i,j)=ci* cj")

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import logm
from user_input import *

###########################################################################
## Solves for eigenfn, eigenenergies of the electronic Hamiltonian at a fixed 
## Calculates time derivative matrix T_der, potential energy: pot_en and acceleration: acc
def tise(x,phi_old):
  Hamil,grad_H=pot(x)
  Ei,phi=np.linalg.eigh(Hamil)

  for i in range(nquant):
    if(sum(phi[:,i]*phi_old[:,i])<0.0):
      phi[:,i]=-phi[:,i]

  T_der=compute_T_der(phi,phi_old)
  phi_old=phi

  pot_en=Ei[lamda]
  for i in range(nclass):
    acc[i]=-1.0/mass[i]*sum(phi[:,lamda]*np.matmul(grad_H[:,:,i],phi[:,lamda]))

  return(pot_en,acc,T_der,phi,phi_old,Ei,grad_H)
###########################################################################

## Classical evolution (x,v) on dtc time step
def evolve_classical(x,v,acc,phi_old):
  ## Velocity Verlet scheme
  ## https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
  x=x+v*dtc+0.5*acc*dtc*dtc
  v=v+0.5*acc*dtc
  pot_en,acc,T_der,phi,phi_old,Ei,grad_H=tise(x,phi_old)
  v=v+0.5*acc*dtc
  energy=pot_en+0.5*sum(mass*v*v)
  return(x,v,acc,energy,phi_old,T_der,Ei,grad_H)
###########################################################################

## Calculates Time derivative matrix T_der as log(U)/dtc
## U=overlap matrix
def compute_T_der(phi,phi_old):
  U=np.zeros((nquant,nquant))
  for i in range(nquant):
    for j in range(nquant):
      U[i,j]=sum(phi_old[:,i]*phi[:,j])
  T_der=logm(U)

  for i in range(nquant):
    for j in range(nquant):
      T_der[i,j]=T_der[i,j]/dtc
  return(T_der)
###########################################################################

## Eqn of d ci/dt
def cidot(ci,Ei,T_der):
  cidot=np.zeros(nquant,dtype=complex)
  for i in range(nquant):
    cidot[i] = Ei[i]*ci[i]/(1.j*hbar)-sum(T_der[i,:]*ci)
  return(cidot)
###########################################################################

## Evolves ci's by a time step dt using RK4 algorithm
## https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
def evolve_rk4(ci,Ei,dt,T_der):
  k1=cidot(ci,Ei,T_der)
  k2=cidot(ci+dt*k1/2,Ei,T_der)
  k3=cidot(ci+dt*k2/2,Ei,T_der)
  k4=cidot(ci+dt*k3,Ei,T_der)
  ci=ci+dt/6*(k1+2*k2+2*k3+k4)
  return(ci)
###########################################################################

## Calculates hopping probability to all states
## if hop_prob[k]<0, hop_prob[k]=0
def compute_hop_prob(dt,ci,T_der):
  hop_prob=np.zeros(nquant)
  for i in range(nquant):
    if(i!=lamda):
      hop_prob[i]=-2*np.real(T_der[i,lamda]*ci[lamda]*np.conj(ci[i]))/abs(ci[lamda]**2) * dt
      if(hop_prob[i]<0.0):
        hop_prob[i]=0.0
  return(hop_prob)
###########################################################################

## Returns state k to which hop is happening
## If hop event does not occur, k is set to 0
def check_hop(hop_prob):
  k=0
  r=np.random.uniform()
  su=0.0
  for i in range(nquant):
    su=su+hop_prob[i]
    if(r<su):
      k=i
      return(k)
  return(k)
###########################################################################

## Evolves ci's over the quantum time step dtq
## At each times step dtq, calculates hopping probability and check if a hop happens
def evolve_quantum(ci,Ei,T_der):
  nsteps = int(dtc/dtq)
  k=0
  for i in range(nsteps):
    ci=evolve_rk4(ci,Ei,dtq,T_der)
    if(k==0):
      hop_prob=compute_hop_prob(dtq,ci,T_der)
      k=check_hop(hop_prob)   ## Returns k=0 if hop does not occur.
  return(k,ci,hop_prob)
###########################################################################

## Calculating pop as average of fraction of trajectories on each state
## rho is the density matrix calculated naively as ci*.cj
## CAREFUL: ci*.cj can give wrong long time answers
def compute_averages(k,ci,lamda):
  pop[lamda,k]+=1
  for i in range(nquant):
    for j in range(nquant):
      rho[i,j,k]+=abs(np.conj(ci[i])*ci[j])
  return(pop,rho)
###########################################################################

## Derivative matrix vector: needed only while hopping
def compute_dij(i,j,phi,grad_H,Ei):
  dij=np.zeros(nclass)
  for k in range(nclass):
    dij[k]=sum(phi[:,i]*np.matmul(grad_H[:,:,k],phi[:,j]))/(Ei[j]-Ei[i])
  return(dij)
###########################################################################

## If a hop event is triggered, this adjusts velocity to conserve energy
## IF HOP IS FRUSTRATED, THIS DOES NOT REVERSE VELOCITY. THIS IS NOT IDEALLY RECOMENDED.
def perform_hop(k,lamda,v,Ei):
  dij=compute_dij(lamda,k,phi,grad_H,Ei)
  aa=sum(0.5*dij*dij/mass)
  bb=sum(v*dij)
  cc=Ei[k]-Ei[lamda]
  if(bb**2-4*aa*cc<0.0):
    ## FRUSTRATED HOP
    ## VELOCITY NEVER REVERSED
    ## INSERT YOUR VELOCITY REVERSAL SCHEME HERE
    gama=0.0
    ## gama = bb/aa ## Will always reverse velocity along dij direction on frustrated hop
    return(lamda,v)
  else:
    if(bb>0.0):
      gama=(-bb+np.sqrt(bb**2-4*aa*cc))/(2*aa)
    else:
      gama=(-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa)
  lamda=k
  v=v-gama*dij/mass
  return(lamda,v)
###########################################################################

def decoherence(ci):
  ## INSERT YOUR DECOHERENCE SCHEME HERE
  return(ci)
###########################################################################

## Evolution by one classical time step dtc
def evolve(x,v,acc,lamda,ci,i,phi_old):
  pop,rho=compute_averages(i,ci,lamda)
  x,v,acc,energy,phi_old,T_der,Ei,grad_H=evolve_classical(x,v,acc,phi_old)
  k,ci,hop_prob=evolve_quantum(ci,Ei,T_der)
  if(k>0):
    lamda,v=perform_hop(k,lamda,v,Ei)
    pot_en, acc, T_der, phi, phi_old, Ei, grad_H = tise(x, phi_old)
    energy = pot_en + 0.5 * sum(mass * v * v)

  ci=decoherence(ci)
  ## Prints information for debugging if nprint set to 1 in user_input.py
  if(nprint==1):
    print(i*dtc,energy,sum(abs(ci**2)),lamda,x,abs(ci[0]**2))
  return(pop,rho,x,v,acc,lamda,ci,phi_old)
###########################################################################

## Prints/plots averages after all the trajecties are done
def averages(pop,rho):
  pop=pop/ntraj
  rho=rho/ntraj
  tim=np.zeros(ntim)
  for i in range(ntim):
    tim[i]=i*dtc
    print(tim[i],pop[:,i],file=pop_file)
    print(tim[i],abs(rho[:,:,i]),file=rho_file)
  plt.plot(tim,pop[0,:])
  plt.plot(tim,abs(rho[0,0,:]))
  plt.plot(tim,abs(rho[1,1,:]))
  plt.plot(tim,abs(rho[0,1,:]))
  plt.show()
  return()
###########################################################################

## Defining required variables
acc=np.zeros(nclass)
phi=np.zeros((nquant,nquant))

pop=np.zeros((nquant,ntim))
rho=np.zeros((nquant,nquant,ntim),dtype=complex)

## Opening files
output = open('energy.out', 'w')
pop_file = open('pop.out', 'w')
rho_file = open('rho.out', 'w')

for i in range(ntraj):
  x,v,lamda,ci=init_cond()    ## init_cond() function in user_input.py
  print("traj number,x,v",i,x,v)
  ## Initializing variables at t=0. Specifically, phi_old and acc are necessary for evolution.
  Hamil,grad_H=pot(x)
  Ei,phi_old=np.linalg.eigh(Hamil)
  pot_en,acc,T_der,phi,phi_old,Ei,grad_H = tise(x,phi_old)
  for j in range(ntim):
    pop,rho,x,v,acc,lamda,ci,phi_old=evolve(x,v,acc,lamda,ci,j,phi_old)

averages(pop,rho)
