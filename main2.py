# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:32:52 2016

@author: Pierre-Edouard
"""


import numpy as np
import matplotlib.pyplot as plt
from SP_FP import *


#equation = 'ADVECTION'; # BURGER or ADVECTION
N = 50; # taille du maillage
Niter = 1; # nombre d'itérations
dt = 0.005; # pas de calcul
c =1.; # vitesse du son
p = 2; #degre polynome dans chaque cellule


## déclaration des vecteurs positions
Xs = defSolPoints(p+1);
Xf = defFluxPoints(p+2);

#déclaration des maillages
Ms = np.zeros((N,p+1));
Mf = np.zeros((N,p+2));

for i in range(N):
    Ms[i,:] = Xs + 2*i +1;
    Mf[i,:] = Xf + 2*i +1;
    
    
##Declaration solution 
U = np.zeros((N,p+1));

##INITIALISATION DE LA SOLUTION
#solution constante
#U = np.ones((N,p+1));
#solution "CHAPEAU"
U[round(N/2),:] = Xs + 1;
U[round(N/2)+1,:] = -Xs + 1;

##declaration matrice Lagrange solution (deplacement vers les points flux)
print(lagrange(Xs,Xf,2,3));


#Matrice d'extrapolation
E = extrapolation(Xs,Xf);

#Matrice de calcul de flux
F = np.zeros((p+2,p+4));
for i in range(1,p+1):
    F[i,i+1]=1.;
F[0,0:2]=[0.5*(1+np.sign(c)),0.5*(1-np.sign(c))];
F[p+1,p+2:p+4]=[0.5*(1+np.sign(c)),0.5*(1-np.sign(c))];
F = c*F;

#déclaration de la matrice U_t "U tild"
U_t = np.zeros((3*(p+1),1));

####################
#début des itérations en temps
for timestep in range(Niter):
    #Début des itérations sur les cellules
    for cell in range(N):
        
        #définition de la variable U_t "U tild"
        if(cell == 0):
            U_t[0:p+1,0] = U[N-1,:];
            U_t[p+1:2*(p+1),0] = U[0,:];
            U_t[2*(p+1):3*(p+1),0] = U[1,:];
        elif(cell == N-1):
            U_t[0:p+1,0] = U[N-2,:];
            U_t[p+1:2*(p+1),0] = U[N-1,:];
            U_t[2*(p+1):3*(p+1),0] = U[0,:];
        else:
            U_t[0:p+1,0] = U[cell-1,:];
            U_t[p+1:2*(p+1),0] = U[cell,:];
            U_t[2*(p+1):3*(p+1),0] = U[cell+1,:];
            
        V = np.dot(E,U_t);
        W = np.dot(F,V);
        
        ###########
        
        ###########
        U[cell,:] = U_t[p+1:2*(p+1),0]        
        
    plt.plot(np.reshape(Ms,(N*(p+1),1)),np.reshape(U,(N*(p+1),1)),'-b')
    plt.show()
    plt.pause(0.0001)