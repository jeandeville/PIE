# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 22:24:19 2015

@author: Etienne
"""




#########################################
#   EQUATION D'ADVECTION DE LA GRANDEUR P
#   dP/dt + c*dP/dx = 0
#   EQUATION DE BURGER DE LA GRANDEUR P
#   dP/dt + c*d(P²)/dx = 0
#########################################


import numpy as np
import matplotlib.pyplot as plt
import fonction as fct

equation = 'ADVECTION'; # BURGER or ADVECTION
N = 10; # taille du maillage
Niter = 200; # nombre d'itérations
dt = 0.005; # pas de calcul
c =1.; # vitesse du son

## déclaration des vecteurs donnés
Pprevious = np.zeros((N,3)); # champ de pression (points solutions) au temps t-dt
P = np.zeros((N,3)); # champ de pression (points solutions)


## déclaration des vecteurs positions
Xs = np.zeros((N,3));
Xf = np.zeros((N,4));

#degre RK
deg = 2

#stockage Qbarre
#Qb1 = np.zeros(deg,N);
#Qb2 = np.zeros(deg,N);
#Qb3 = np.zeros(deg,N);


## initialisation des vecteurs position
for i in range(N):
    for j in range(3):
        Xs[i,j] = (2.*j+1.)/6./N + i*1./N;
        
for i in range(N):
    for j in range(4):
        Xf[i,j] = j/3./N+ i*1./N;

## initialisation du vecteur de pression
for cell in range(np.round(N/2)-1, np.round(N/2)):#############################CONDITION INITIALES
    for j in range(3):
        P[cell,j] = np.sin(np.pi*N*Xs[cell,j]);



for iter in range(Niter):
    
    Pprevious = P;
    
        #RK1 
    P = Pprevious - dt * fct.flux(Pprevious,Xs,Xf,N,c)
    
    
        #RK2
    #for i in range(deg):
        
      

    
        # plot du champ P
    plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(P,(3*N,1)),'-b')
    plt.show()
    plt.pause(0.0001)

    
    





