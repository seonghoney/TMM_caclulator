# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:01:18 2024

@author: Seongheon Kim

(input)
n : array of refractive indices for each layer
d : array of thicknesses of each layer
w : frequency

"""
import numpy as np
from numpy.linalg import inv
from math import pi

c = 299792458 # Speed of Light

def T_matrix_calculator(n, d, w, theta = 0):
    
    k1 = n*w/c*np.cos(theta/180*pi)
    T_matrix = np.zeros([len(w),2,2],dtype=np.complex_)
    
    for i in range(len(w)):
        T_matrix[i][:][:] = np.matrix([[np.cos(k1[i]*d), -1j*np.sin(k1[i]*d)/n],
                                       [-1j*np.sin(k1[i]*d)*n, np.cos(k1[i]*d)]])
       
    return T_matrix

def solve(matrice, profile):
    for i in range(len(profile)):
        if i == 0:
            TMM = matrice[profile[0]]
        else:
            TMM = matrice[profile[i]]@TMM # [# of freqs, 2, 2]
            
    bandstrure_k = np.arccos(0.5*np.trace(TMM,offset=0, axis1=1, axis2=2))
   
    return TMM, bandstrure_k

def RT_coeff(w, TMM, N1):
    
    transmission = np.zeros([len(w)], dtype = np.complex_)
    reflection = np.zeros([len(w)], dtype = np.complex_)
    
    for i in range(len(w)):
        cal_M = inv(N1)@inv(TMM[i][:][:])@N1
        transmission[i] = 1/cal_M[0,0]
        reflection[i] = cal_M[1,0]*transmission[i]
        
    return reflection, transmission 
