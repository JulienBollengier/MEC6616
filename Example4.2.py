# -*- coding: utf-8 -*-
""" LAP1 """

""" Example 4.1 """
""" Source heat conduction in an insulated rod """
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt

""" Constants """
L = 0.02; # Rod distance in meters [m] 
n = 5; # Number of nodes [nb]
k = 0.5; # Thermal conductivity [W/m.K]
q = 1000e3; # Heat Generation [W/m^3]
A = 1; # Cross-section Area [m^2]
Ta = 373.15; # Temperature at the beginning in degres [K]
Tb = 473.15; # Temperature at the end in degres [K]

""" Calculated constants """
dx = L/n; # Distance between nodes [m] """
TaD = Ta - 273.15; # Temperature at the beginning in degres [°C] 
TbD = Tb - 273.15; # Temperature at the end in degres [°C] 


""" Discretised equation for nodal points 2,3 & 4 : apTp = awTw + aeTe + Su"""
aw = k*A/dx;
ae = k*A/dx;
Sp = 0;
Su = q*A*dx;
ap = aw + ae - Sp;

""" Discretised equation for nodal point 1 (boundary nodes) : apTp = awTw + aeTe + Su"""
aw1 = 0;
ae1 = k*A/dx;
Su1 = q*A*dx + 2*k*A/dx*Ta;
Sp1 = -2*k*A/dx
ap1 = aw1 + ae1 - Sp1;

""" Discretised equation for nodal point 5 (boundary nodes) : apTp = awTw + aeTe + Su"""
aw5 = k*A/dx;
ae5 = 0;
Su5 = q*A*dx + 2*k*A/dx*Tb;
Sp5 = -2*k*A/dx
ap5 = aw5 + ae5 - Sp5;



""" Equation's system assembly """
i = 0;
j = dx/2;
X = np.zeros((n+2,1));
M = np.zeros((n,n));
S = np.zeros((n,1));

while i < n :
    if i == 0 :
        M[i,i] = ap1
        M[i,i+1] = -ae1
        S[i,:] = Su1
        X[i+1,:] = j
    if i == n-1 :
        M[i,i] = ap5
        M[i,i-1] = -aw5
        S[i,:] = Su5
        X[i+2,:] = j + dx/2
    if i > 0 and i < n-1:
        M[i,i] = ap
        M[i,i+1] = -aw
        M[i,i-1] = -ae
        S[i,:] = Su
        X[i+1,:] = j
        X[i+2,:] = j + dx
    i = i+1
    j = j + dx

""" System Solutions """
T = solve(M,S)-273.15;


""" Graphs """
""" Exact Solution : Ttheo = [(Tb-Ta)/L + q/2k(L-x)]x + Ta """
X = np.linspace(0,L,n+2);
Y = [[TaD]] + T.tolist() + [[TbD]];
x = np.linspace(0,L,50);
y = ((TbD-TaD)/L + q/(2*k)*(L-x))*x + TaD;


plt.grid()
plt.title('Comparaison of the Num. res. w/ the Analytical sol.');
plt.xlabel('Distance [m]');
plt.ylabel('Temperature [°C]');
plt.plot(X,Y,'r+',x,y,'b--');
plt.legend(['Numerical','Analytical']);
plt.show()

