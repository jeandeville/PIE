# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:25:32 2016

@author: jeandeville
"""

import numpy as np

def jacobian(Xirreg):
    
    N = np.size(Xirreg) - 1;
    
    jac = np.zeros([N,3]);
    
    for i in range(N):
        jac[i,0] = (Xirreg[i+1] - Xirreg[i])/2;
        jac[i,1] = (Xirreg[i+1] + Xirreg[i])/2;
        jac[i,2] = 1/jac[i,0];
        
    return jac;



def maillage(Xs, jac,N,p):
    

    
    Ms_irreg = np.zeros([N,p+1]);
    
    for i in range(N):
        Ms_irreg[i,:] = np.add(np.multiply(jac[i,0],Xs), jac[i,1]); 
        
    return Ms_irreg;