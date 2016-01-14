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


def flux(Q,Xs,Xf,N,c):
    
    equation ='ADVECTION'

    #N taille du maillage
    #c vitesse du son

    ## déclaration des vecteurs donnés
    
    F = np.zeros((N,4)); # champ de fux (points flux)
    G = np.zeros((N,3)); # G = div(F)
            
    ## déclaration des vecteurs de Lagrange (utilisés pour stocker les coefs ...
    #        des polynomes portant le même nom)
    Ls = np.zeros((N,3));        
    Lf = np.zeros((N,4));
    Lg = np.zeros((N,3));

    
    for cell in range(N):
        
        # calcul des coef du polynome de Lagrange pour P
        a1 = 0.;
        b1 = 0.;
        c1 = 0.;
        for i in range (3):
            sump = 0.;
            prodp = 1.;
            denom = 1.;
            for j in range(3):
                if (i != j):
                    sump = sump + Xs[cell,j];
                    prodp = prodp*Xs[cell,j];
                    denom = denom*(Xs[cell,i]-Xs[cell,j]);
            a1 = a1 + Q[cell,i]/denom;
            b1 = b1 - Q[cell,i]*sump/denom;
            c1 = c1 + Q[cell,i]*prodp/denom;
        Ls[cell,0] = a1;
        Ls[cell,1] = b1;
        Ls[cell,2] = c1;
        
        # interpolation de P aux points flux
        for j in range(4):
            F[cell,j] = Ls[cell,0]*Xf[cell,j]**2.+Ls[cell,1]*Xf[cell,j]+Ls[cell,2];

            # reconstruire la véritable valeur de F (euler -> non linéaire)
            if equation == 'ADVECTION':
                F[cell,j] = c*F[cell,j];
            if equation == 'BURGER':
                F[cell,j] = c*F[cell,j]*F[cell,j];
                
    # écriture du solveur de Riemann (différentiation en fonction du signe de c)
    if c >= 0:
        F[0,0] = F[N-1,3];
        for cell in range(1,N):
            F[cell,0] = F[cell-1,3]
    else:
        F[N-1,3] = F[0,0];
        for cell in range(1,N):
            F[cell-1,3] = F[cell,0]
              
    # calcul de la contribution des points flux sur les points solution
    for cell in range(N):
        # calcul des coef du polynome de Lagrange pour F
        a2 = 0.;
        b2 = 0.;
        c2 = 0.;
        d2 = 0.;
        for i in range (4):
            sump = 0.;
            prodp = 1.;
            prodp3 = 0.;
            denom = 1.;
            for j in range(4):
                if (i != j):
                    sump = sump + Xf[cell,j];
                    prodp = prodp*Xf[cell,j];
                    denom = denom*(Xf[cell,i]-Xf[cell,j]);
                # calcul de la somme des doubles produits (= coef du terme en "x")
                    prodp2 = 1;             
                    for k in range(4):
                        if (k != i)&(k != j):
                            prodp2 = prodp2*Xf[cell,k];
                    prodp3 = prodp3 + prodp2;
                
            a2 = a2 + F[cell,i]/denom;
            b2 = b2 - F[cell,i]*sump/denom;
            c2 = c2 + F[cell,i]*prodp3/denom;
            d2 = d2 - F[cell,i]*prodp/denom; 
        Lf[cell,0] = a2;
        Lf[cell,1] = b2;
        Lf[cell,2] = c2;
        Lf[cell,3] = d2;
        
        # calcul de l'interpolation de la divergence du flux G = dF/dx
        Lg[cell,0] = 3.*Lf[cell,0];
        Lg[cell,1] = 2.*Lf[cell,1];
        Lg[cell,2] = Lf[cell,2];
        
        # évaluation de G aux points solution
        for j in range(3):
            G[cell,j] = Lg[cell,0]*Xs[cell,j]**2.+Lg[cell,1]*Xs[cell,j]+Lg[cell,2];

    return G
