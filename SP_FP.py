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

def adaptSolPoints(N):
    Xp = np.zeros((N,3));
    xxp=defSolPoints(N*(3-1)+1);
    #print("xxp",xxp);
    
    for i in range(N):
            
            #print("Ligne", i, Xp[i,:])
            #print("Ligne", i, xxp[i*2:(i+1)*2+1])
            Xp[i,:] = xxp[i*2:(i+1)*2+1];
         
    #print("Xp",Xp);
    return Xp;
    
def adaptFluxPoints(N):
    Xf = np.zeros((N,4));
    #print("Xf",Xf);
    xxf=defFluxPoints(N*(4-1)+1);
    #print("xxf",xxf);
    
    for i in range(N):
            
            #print("Ligne", i, Xf[i,:])
            #print("Ligne", i, xxf[i*3:(i+1)*3+1])
            Xf[i,:] = xxf[i*3:(i+1)*3+1];
         
    #print("Xf",Xf);
    return Xf;
    
    
    
    