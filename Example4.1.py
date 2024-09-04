# -*- coding: utf-8 -*-
""" LAP1 """

""" Example 4.1 """
""" Source-free heat conduction in an insulated rod """
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt

""" Constants """
L = 0.5; # Rod distance in meters [m]
n = 5; # Number of nodes [nb] """
k = 1000; # Thermal conductivity [W/m.K]
A = 10e-3; # Cross-section Area [m^2]
Ta = 100+273.15; # Temperature at the beginning in degres [K]
Tb = 500+273.15; # Temperature at the end in degres [K]

""" Calculated constants """
dx = L/n; # Distance between nodes [m]
Sp = (-2)*k*A/dx;
Su = 2*k*A/dx*0;


""" Discretised equation for nodal points 2,3 & 4 : apTp = awTw + aeTe"""
aw = k*A/dx;
ae = k*A/dx;
ap = aw + ae; 

""" Discretised equation for nodal point 1 (boundary nodes) : apTp = awTw + aeTe + Su"""
aw1 = 0;
ae1 = k*A/dx;
ap1 = aw1 + ae1 - Sp;
Su1 = 2*k*A/dx*Ta;

""" Discretised equation for nodal point 5 (boundary nodes) : apTp = awTw + aeTe + Su"""
aw5 = k*A/dx;
ae5 = 0;
ap5 = aw5 + ae5 - Sp;
Su5 = 2*k*A/dx*Tb;



""" Equation's system assembly """ 
i = 0;
M = np.zeros((n,n));
S = np.zeros((n,1));

while i < n :
    if i == 0 :
        M[i,i] = ap1
        M[i,i+1] = -ae1
        S[i,:] = Su1
    if i == n-1 :
        M[i,i] = ap5
        M[i,i-1] = -aw5
        S[i,:] = Su5
    if i > 0 and i < n-1:
        M[i,i] = ap
        M[i,i+1] = -aw
        M[i,i-1] = -ae
        S[i,:] = Su
    i=i+1


""" System Solutions """
T = solve(M,S)-273.15;


""" Graphs """
""" Exact Solution : Ttheo = 800x + 100 """
X = np.linspace(0,L,n+2);
Y = [[Ta-273.15]] + T.tolist() + [[Tb-273.15]];
x = np.linspace(0,L,50);
y = 800*x + 100;

""" Errors """


""" Graphics Display """
plt.grid()
plt.title('Comparaison of the Num. res. w/ the Analytical sol.');
plt.xlabel('Distance [m]');
plt.ylabel('Temperature [°C]');
plt.plot(X,Y,'r+',x,y,'b--');
plt.legend(['Numerical','Analytical']);
plt.show()

#test
