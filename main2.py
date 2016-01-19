# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:32:52 2016
@author: jeandeville
"""


import numpy as np
import matplotlib.pyplot as plt
from SP_FP import *


#equation = 'ADVECTION'; # BURGER or ADVECTION
N = 50; # taille du maillage
Niter = 500; # nombre d'itérations
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
U = np.zeros(N*(p+1));

##declaration matrice Lagrange solution (deplacement vers les points flux)
#print(lagrange(Xs,Xf,2,3));

##Matrice 
E = extrapolation(Xs,Xf);
print("matrice d'éxtrapolation",E);