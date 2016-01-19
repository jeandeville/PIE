# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 19:55:54 2016

@author: etienne
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 22:24:19 2015

@author: Etienne
"""




#####################################################
#                EQUATION D'EULER 1D
#
#                 dQ/dt + dF/dx = 0
#
#       Q = [rho, rho*u, rho*E]^T
#       F = [rho*u, P + rho*u^2, u*(rho*E+P) ]^T
#            
#            [                 Q_2                         ]
#       F =  [     Kp*Q_3 + (1-0.5*Kp)*(Q_2^2/Q_1)         ]
#            [ (1+Kp)*(Q_2*Q_3/Q_1)-0.5*Kp*(Q_2^3/Q_1^2)   ]
#
#                   Kp = gama -1
#####################################################


import numpy as np
import matplotlib.pyplot as plt

equation = 'ADVECTION'; # BURGER or ADVECTION
N = 10; # taille du maillage
Niter = 100; # nombre d'itérations
dt = 0.001; # pas de calcul

gama = 1.4; # constante du GP
Kp = gama-1;

c = 1.; # sens des conditions de Riemann

## déclaration des vecteurs donnés
Q_1 = np.zeros((N,3)); # champ de pression (points solutions)
Q_2 = np.zeros((N,3)); #
Q_3 = np.zeros((N,3)); # 
F_1 = np.zeros((N,4)); # champ de fux (points flux)
F_2 = np.zeros((N,4)); #
F_3 = np.zeros((N,4)); #
F_1storage = np.zeros((N,4)); # champ de fux (points flux)
F_2storage = np.zeros((N,4)); #
F_3storage = np.zeros((N,4)); #
G_1 = np.zeros((N,3)); # G = div(F)
G_2 = np.zeros((N,3)); # 
G_3 = np.zeros((N,3)); # 
H = np.zeros((N,3)); # vecteur test pour vérifier l'interpolation Lagr4


## déclaration des vecteurs positions
Xs = np.zeros((N,3));
Xf = np.zeros((N,4));

## initialisation des vecteurs position
for i in range(N):
    for j in range(3):
        Xs[i,j] = (2.*j+1.)/6./N + i*1./N;
        
for i in range(N):
    for j in range(4):
        Xf[i,j] = j/3./N+ i*1./N;

## initialisation du vecteur de pression

for cell in range(0, N):#############################CONDITION INITIALES
    for j in range(3):
        Q_1[cell,j] = 1.;
        Q_3[cell,j] = 0.5;
for cell in range(np.round(N/2)-1, np.round(N/2)):
    for j in range(3):
        Q_2[cell,j] = np.sin(np.pi*N*Xs[cell,j]);

Q_1previous = Q_1;
Q_2previous = Q_2;
Q_3previous = Q_3;
        
## déclaration des vecteurs de Lagrange (utilisés pour stocker les coefs ...
#        des polynomes portant le même nom)
Lq_1 = np.zeros((N,3));     
Lq_2 = np.zeros((N,3));  
Lq_3 = np.zeros((N,3));

Lf_1 = np.zeros((N,4));
Lf_2 = np.zeros((N,4));
Lf_3 = np.zeros((N,4));

Lg_1 = np.zeros((N,3));
Lg_2 = np.zeros((N,3));
Lg_3 = np.zeros((N,3));

## début des itérations
for t in range(Niter):
    
    for cell in range(N):
        
        # calcul des coef du polynome de Lagrange pour P
        a1 = np.zeros(3);
        b1 = np.zeros(3);
        c1 = np.zeros(3);
        for i in range (3):
            sump = 0.;
            prodp = 1.;
            denom = 1.;
            for j in range(3):
                if (i != j):
                    sump = sump + Xs[cell,j];
                    prodp = prodp*Xs[cell,j];
                    denom = denom*(Xs[cell,i]-Xs[cell,j]);
                    
            a1[0] = a1[0] + Q_1[cell,i]/denom;
            a1[1] = a1[1] + Q_2[cell,i]/denom;
            a1[2] = a1[2] + Q_3[cell,i]/denom;
            
            b1[0] = b1[0] - Q_1[cell,i]*sump/denom;
            b1[1] = b1[1] - Q_2[cell,i]*sump/denom;
            b1[2] = b1[2] - Q_3[cell,i]*sump/denom;
            
            c1[0] = c1[0] + Q_1[cell,i]*prodp/denom;
            c1[1] = c1[1] + Q_2[cell,i]*prodp/denom;
            c1[2] = c1[2] + Q_3[cell,i]*prodp/denom;
            
        Lq_1[cell,0] = a1[0];
        Lq_1[cell,1] = b1[0];
        Lq_1[cell,2] = c1[0];
        
        Lq_2[cell,0] = a1[1];
        Lq_2[cell,1] = b1[1];
        Lq_2[cell,2] = c1[1];
        
        Lq_3[cell,0] = a1[2];
        Lq_3[cell,1] = b1[2];
        Lq_3[cell,2] = c1[2];
        
        # interpolation de P aux points flux
        for j in range(4):
            F_1storage[cell,j] = Lq_1[cell,0]*Xf[cell,j]**2.+Lq_1[cell,1]*Xf[cell,j]+Lq_1[cell,2];
            F_2storage[cell,j] = Lq_2[cell,0]*Xf[cell,j]**2.+Lq_2[cell,1]*Xf[cell,j]+Lq_2[cell,2];
            F_3storage[cell,j] = Lq_3[cell,0]*Xf[cell,j]**2.+Lq_3[cell,1]*Xf[cell,j]+Lq_3[cell,2];

            # reconstruire la véritable valeur de F (euler -> non linéaire)           
            #            [                 Q_2                         ]
            #       F =  [     Kp*Q_3 + (1-0.5*Kp)*(Q_2^2/Q_1)         ]
            #            [ (1+Kp)*(Q_2*Q_3/Q_1)-0.5*Kp*(Q_2^3/Q_1^2)   ]
            
            F_1[cell,j] = F_2storage[cell,j];
            F_2[cell,j] = Kp*F_3storage[cell,j]+ (1-0.5*Kp)*(F_2storage[cell,j]**2/F_1storage[cell,j]);
            F_3[cell,j] = (1+Kp)*(F_2storage[cell,j]*F_3storage[cell,j]/F_1storage[cell,j])-0.5*Kp*(F_2storage[cell,j]**3/F_1storage[cell,j]**2)
                
    # écriture du solveur de Riemann (différentiation en fonction du sugne de c)
    if c >= 0:
        F_1[0,0] = F_1[N-1,3];
        F_2[0,0] = F_2[N-1,3];
        F_3[0,0] = F_3[N-1,3];
        
        for cell in range(1,N):
            F_1[cell,0] = F_1[cell-1,3];
            F_2[cell,0] = F_2[cell-1,3];
            F_3[cell,0] = F_3[cell-1,3];
    else:
        F_1[N-1,3] = F_1[0,0];
        F_2[N-1,3] = F_2[0,0];
        F_3[N-1,3] = F_3[0,0];
        for cell in range(1,N):
            F_1[cell-1,3] = F_1[cell,0];
            F_2[cell-1,3] = F_2[cell,0];
            F_3[cell-1,3] = F_3[cell,0];
              
    # calcul de la contribution des points flux sur les points solution
    for cell in range(N):
        # calcul des coef du polynome de Lagrange pour F
        a2 = np.zeros(3);
        b2 = np.zeros(3);
        c2 = np.zeros(3);
        d2 = np.zeros(3);
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
                
            a2[0] = a2[0] + F_1[cell,i]/denom;
            a2[1] = a2[1] + F_2[cell,i]/denom;
            a2[2] = a2[2] + F_3[cell,i]/denom;
            
            b2[0] = b2[0] - F_1[cell,i]*sump/denom;
            b2[1] = b2[1] - F_2[cell,i]*sump/denom;
            b2[2] = b2[2] - F_3[cell,i]*sump/denom;
            
            c2[0] = c2[0] + F_1[cell,i]*prodp3/denom;
            c2[1] = c2[1] + F_2[cell,i]*prodp3/denom;
            c2[2] = c2[2] + F_3[cell,i]*prodp3/denom;
                # permet de tester l'interpolation Lagr3
    #plt.plot(np.reshape(Xf,(4*N,1)),np.reshape(F_test,(4*N,1)),'-r')
    #plt.plot(np.reshape(Xp,(3*N,1)),np.reshape(Pprevious,(3*N,1)),'.')
            d2[0] = d2[0] - F_1[cell,i]*prodp/denom;
            d2[1] = d2[1] - F_2[cell,i]*prodp/denom; 
            d2[2] = d2[2] - F_3[cell,i]*prodp/denom; 
            
        Lf_1[cell,0] = a2[0];
        Lf_1[cell,1] = b2[0];
        Lf_1[cell,2] = c2[0];
        Lf_1[cell,3] = d2[0];
        
        Lf_2[cell,0] = a2[1];
        Lf_2[cell,1] = b2[1];
        Lf_2[cell,2] = c2[1];
        Lf_2[cell,3] = d2[1];

        Lf_3[cell,0] = a2[2];
        Lf_3[cell,1] = b2[2];
        Lf_3[cell,2] = c2[2];
        Lf_3[cell,3] = d2[2];
        
        for j in range(3):
            H[cell,j] = Lf_2[cell,0]*Xs[cell,j]**3.+Lf_2[cell,1]*Xs[cell,j]**2.+Lf_2[cell,2]*Xs[cell,j]+Lf_2[cell,3];
        
        # calcul de l'interpolation de la divergence du flux G = dF/dx
        Lg_1[cell,0] = 3.*Lf_1[cell,0];
        Lg_1[cell,1] = 2.*Lf_1[cell,1];
        Lg_1[cell,2] = Lf_1[cell,2];
        
        Lg_2[cell,0] = 3.*Lf_2[cell,0];
        Lg_2[cell,1] = 2.*Lf_2[cell,1];
        Lg_2[cell,2] = Lf_2[cell,2];
        
        Lg_3[cell,0] = 3.*Lf_3[cell,0];
        Lg_3[cell,1] = 2.*Lf_3[cell,1];
        Lg_3[cell,2] = Lf_3[cell,2];
        
        # évaluation de G aux points solution
        for j in range(3):
            G_1[cell,j] = Lg_1[cell,0]*Xs[cell,j]**2.+Lg_1[cell,1]*Xs[cell,j]+Lg_1[cell,2];
            G_2[cell,j] = Lg_2[cell,0]*Xs[cell,j]**2.+Lg_2[cell,1]*Xs[cell,j]+Lg_2[cell,2];
            G_3[cell,j] = Lg_3[cell,0]*Xs[cell,j]**2.+Lg_3[cell,1]*Xs[cell,j]+Lg_3[cell,2];
            
    # avancement d'un pas de temps
    Q_1storage = Q_1;
    Q_2storage = Q_2;
    Q_3storage = Q_3;
    
    #Bete avancement en Euler
    Q_1 = Q_1 - G_1*dt;
    Q_2 = Q_2 - G_2*dt;
    Q_3 = Q_3 - G_3*dt;
    Q_1previous = Q_1storage;    
    Q_2previous = Q_2storage;  
    Q_3previous = Q_3storage;  
    
    # plot des champs
    plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(Q_1,(3*N,1)),'-b')
    plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(Q_2,(3*N,1)),'-r')
    plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(Q_3,(3*N,1)),'-g')
    
    # permet de tester l'interpolation Lagr3
    #plt.plot(np.reshape(Xf,(4*N,1)),np.reshape(F_2storage,(4*N,1)),'-r')
    #plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(Q_2,(3*N,1)),'.')
    
    # permet de tester l'interpolation Lagr4
    #plt.plot(np.reshape(Xs,(3*N,1)),np.reshape(H,(3*N,1)),'-r')
    #plt.plot(np.reshape(Xf,(4*N,1)),np.reshape(F_2,(4*N,1)),'.')
    
    
    plt.show()
    plt.pause(0.01)
    
print(G_1)




