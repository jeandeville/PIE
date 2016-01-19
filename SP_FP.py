#Import modules
import numpy as np
import numpy.polynomial.legendre as leg
import sympy as sp
import matplotlib.pyplot as plt
#Module : mathematical fonctions
import math as ma

#Definitions of solutions points :
def defSolPoints(n):
    res=np.zeros(n);
    for i in range(n):
        res[i]=-ma.cos(((2.0*(i+1)-1)/(2*n))*ma.pi);
    return res;
    
#Definitions of flux points :
def defFluxPoints(n):
    res=np.zeros(n);
    res[0]=-1.0;
    res[-1]=1.0;
    coeff=np.zeros(n-1);
    coeff[-1]=1.0;
    roots=leg.legroots(coeff);
    for i in range(1,n-1):
        res[i]=roots[i-1];
    return res;
    
      
        
# Definition du terme generique de la matrice  [L_k(f_j)](1<=k<=p+1;1<=j<=p+2)
def lagrange(Xs,k,x):
    
    p = np.size(Xs)-1;
    
    Li = 1
    for i in range(p+1):
        if (i!=k):
            Li = Li*(x-Xs[i])/(Xs[k]-Xs[i])
                        
    return Li;
    
    
 # Definition du terme generique de la matrice [T_k(s_j)](1<=j<=p+1)   
def lagrangeTiPrime(Xf,k,x):
    
    p = np.size(Xf)-2;
    Tiprime=0;
    
    for i in range(p+2):
        if (i!=k):
            Tiprime = Tiprime + (1/(Xf[k]-Xf[i]))*lagrange(Xf,k,x)/( (x-Xf[i])/(Xf[k]-Xf[i]) );
            
    return Tiprime;
            
            
  # Creation de la matrice d'extrapolation
def extrapolation(Xs,Xf):
    
    p = np.size(Xs)-1;
    
    E = np.zeros([p+4,3*(p+1)]);
    
    for i in range(p+1):
        E[0,i] = lagrange(Xs,i,Xf[p+1]);
    
    for i in range(p+2):
        for j in range(p+1):
            E[i+1,p+1+j] = lagrange(Xs,j,Xf[i]);
            
    for i in range(p+1):
        E[p+3,2*(p+1)+i] = lagrange(Xs,i,Xf[0]);
    return E;
    
# Creation de la matrice de derivation
def derivation(Xs,Xf):
    
    p = np.size(Xs)-1;
    D = np.zeros([p+1,p+2]);
    
    for i in range(p+1):
        for j in range(p+2):
            D[i,j] = lagrangeTiPrime(Xf,j,Xs[i]);
            
    return D;
    


    