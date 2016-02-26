# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:32:52 2016

@author: Pierre-Edouard
"""


import numpy as np
import matplotlib.pyplot as plt
from SP_FP import *
from fonctions_irreg import *

c =1.; # vitesse du son
CFL = 1.;
p = 3; #degre polynome dans chaque cellule
#equation = 'ADVECTION'; # BURGER or ADVECTION
N = 25; # taille du maillage, nombre de cellules
Niter = 0; # nombre d'itérations
dt = CFL*2/(c*(p+1)); # pas de calcul
###########
###
L = 50.;
####
###########
####################
####################
####################
##########

#Definition

#Xirreg = np.linspace(0,L,N+1);
Xirreg = [0.,3.,4.,6.,9.,11.,12.,13.,15.,19.,20.,24.,25.,26.,28.,32.,33.,34.,38.,42.,43.,45.,46.,47.,48.,50.]; #maillage des cellules

jacobian = jacobian(Xirreg);

Ms_irreg = np.zeros([N,p+1]); #maillage irrégulier des points solutions

##########
####################
####################
####################


## déclaration des vecteurs positions
Xs = defSolPoints(p+1);
Xf = defFluxPoints(p+2);

#déclaration des maillages
Ms = np.zeros((N,p+1));
Mf = np.zeros((N,p+2));

for i in range(N):
    Ms[i,:] = Xs + 2*i +1;
    Mf[i,:] = Xf + 2*i +1;

####################
####################
####################
#########

Ms_irreg = maillage(Xs,jacobian,N,p); #définition de l'emplacement des points solutions dans le "physical plane" à partir de ceux du "computational plane"


U = np.exp(-(np.add(Ms_irreg,-L/2)**2/100)); #initialisation dans le domaine irregulier

#plt.plot(np.reshape(Ms_irreg,(N*(p+1),1)),np.reshape(U,(N*(p+1),1)),'-b')
#plt.show()
    
Ubarre = np.zeros([N,p+1]);

for i in range(p+1):
    Ubarre[:,i] = np.multiply(jacobian[:,0],U[:,i]);
    
plt.plot(np.reshape(Ms,(N*(p+1),1)),np.reshape(Ubarre,(N*(p+1),1)),'-b')
plt.show()

#########
####################
####################
####################   

 
#U = np.exp(-(np.add(Ms,-100)**2/0.1));
#U = np.exp(-((np.add(Ms,-100)+0.)**2/100.));







##Declaration solution 
#U = np.zeros((N,p+1));

##INITIALISATION DE LA SOLUTION
#solution constant
#U = np.ones((N,p+1));
#solution "CHAPEAU"
#U[round(N/2),:] = Xs + 1;
#U[round(N/2)+1,:] = -Xs + 1;



#Matrice d'extrapolation
E = extrapolation(Xs,Xf);

#Matrice de calcul de flux
F = np.zeros((p+2,p+4));
for i in range(1,p+1):
    F[i,i+1]=1.;
F[0,0:2]=[0.5*(1+np.sign(c)),0.5*(1-np.sign(c))];
F[p+1,p+2:p+4]=[0.5*(1+np.sign(c)),0.5*(1-np.sign(c))];
F = c*F;
#creation de la matrice de dérivation
D = derivation(Xs,Xf);
#déclaration de la matrice U_t "U tild"
U_t = np.zeros((3*(p+1),1));

#declaration de M matrice N x (p+1)
M = np.zeros((N,p+1));

#stockage coefficients RK (cf article Bogey 2013)
gamma = np.zeros(6);
gamma[5]=1;
gamma[4]=0.5;
gamma[3]=0.165919771368;
gamma[2]=0.040919732041;
gamma[1]=0.007555704391;
gamma[0]=0.000891421261;

####################
#début des itérations en temps
for timestep in range(Niter):
    
    Ustockage = Ubarre;
    
    #boucle de RK6    
    for i in range(6):    
        
        #Début des itérations sur les cellules
        for cell in range(N):
            
            #définition de la variable U_t "U tild"
            if(cell == 0):
                U_t[0:p+1,0] = Ubarre[N-1,:];
                U_t[p+1:2*(p+1),0] = Ubarre[0,:];
                U_t[2*(p+1):3*(p+1),0] = Ubarre[1,:];
            elif(cell == N-1):
                U_t[0:p+1,0] = Ubarre[N-2,:];
                U_t[p+1:2*(p+1),0] = Ubarre[N-1,:];
                U_t[2*(p+1):3*(p+1),0] = Ubarre[0,:];
            else:
                U_t[0:p+1,0] = Ubarre[cell-1,:];
                U_t[p+1:2*(p+1),0] = Ubarre[cell,:];
                U_t[2*(p+1):3*(p+1),0] = Ubarre[cell+1,:];
                
            V = np.dot(E,U_t);
            W = np.dot(F,V);
            M_t = np.dot(D,W);
            
            ###########      
            M[cell,:] =np.transpose(M_t);     
       
        Mbarre = np.zeros([N,p+1]);
        for j in range(p+1):
            for k in range(N):
                Mbarre[k,j] = M[k,j];#*jacobian[k,2];
            
        Ubarre = Ustockage - gamma[i] * dt * Mbarre;
        for j in range(p+1):
            for k in range(N):
                U[k,j] = Ubarre[k,j]*jacobian[k,2];
        
    
        
        
    plt.plot(np.reshape(Ms_irreg,(N*(p+1),1)),np.reshape(U,(N*(p+1),1)),'-b')
    plt.show()
    plt.pause(0.0000000001)
    