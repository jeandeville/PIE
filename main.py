# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 22:24:19 2015

@author: Etienne

Par rapport à main: suppression de RK1, implémentation de RK2
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
P = np.zeros((N,3)); # champ de pression (points solutions)


## déclaration des vecteurs positions
Xs = np.zeros((N,3));
Xf = np.zeros((N,4));

#degre RK
deg = 2;

#stockage Qbarre
Qb = np.zeros((N,3,deg));



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
        Qb[cell,j,0] = np.sin(np.pi*N*Xs[cell,j]);



for iter in range(Niter):
    
        # plot du champ P
    plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(Qb[:,:,0],(3*N,1)),'-b')
    plt.show()
    plt.pause(0.0001)    
    
        #RK2: calcul des Qbarres
    for i in range(1,deg):
        Qb[:,:,i] = Qb[:,:,0] - 0.5 * dt * fct.flux(P,Xs,Xf,N,c);
        
        #Assemblage final Runge-Kutta
    Qb[:,:,0] = Qb[:,:,0] - dt * fct.flux(Qb[:,:,1],Xs,Xf,N,c);
        
      



    
    





