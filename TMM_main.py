# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:54:00 2024

@author: Seongheon Kim
"""
import numpy as np


import matplotlib.pyplot as plt
from TMM_solver import T_matrix_calculator, solve, RT_coeff


from math import pi

c = 299792458


theta = 45 #incident angle

w = 2*pi*c*np.linspace(0,0.4,401)

n1 = 4.6;
d1 = 1/3;

n2 = 1.6;
d2 = 2/3;

N1 =np.matrix([[1,1],[-1,1]]) # air


T1 = T_matrix_calculator(n1,d1,w, theta)
T2 = T_matrix_calculator(n2,d2,w, theta)

matrice = {1: T1, 2:T2}
profile = [1,2,1,2,1] # top to bottom

TMM, bandstructure = solve(matrice, profile)
r, t = RT_coeff(w, TMM, N1)
 
norm_w = w/(2*pi*c)

plt.plot(norm_w,abs(t)**2, norm_w, abs(r)**2)
plt.axis([0,0.4,0,1])
plt.ylabel('Reflection (Transmission)')
plt.xlabel('w(2$\pi$c/a)')
plt.show()

plt.plot(bandstructure/(2*pi), norm_w)
plt.axis([0,0.5,0,0.4])
plt.xlabel('k(2$\pi$/a)')
plt.ylabel('w(2$\pi$c/a)')
plt.show()