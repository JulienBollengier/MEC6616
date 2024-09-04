# -*- coding: utf-8 -*-
""" LAP1 """

""" Example 4.1 """
""" Cooling of a circular fin by means of convective heat transfer along its length. """
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt

""" Constants """
L = 1; # Rod distance in meters [m]
n = 5; # Number of nodes [nb] 
q = 1000e3; # Heat Generation [W/m^3]
n2 = n**2; # n^2 = h.P/k.A [m^(-2)] 
Ta = 293.15; # Ambiant Temperature in degres [K]
Tb = 373.15; # Temperature at the end in degres [K]

""" Calculated constants """
dx = L/n; # Distance between nodes [m] 
TaD = Ta - 273.15; # Temperature at the beginning in degres [°C] 
TbD = Tb - 273.15; # Temperature at the end in degres [°C] 


""" Discretised equation for nodal points 2,3 & 4 : apTp = awTw + aeTe + Su"""
aw = 1/dx;
ae = 1/dx;
Sp = -n2*dx;
Su = n2*dx*Ta;
ap = aw + ae - Sp;

""" Discretised equation for nodal point 1 (boundary nodes) : apTp = awTw + aeTe + Su"""
aw1 = 0;
ae1 = 1/dx;
Su1 = n2*dx*Ta + 2*Tb/dx;
Sp1 = -n2*dx-2/dx
ap1 = aw1 + ae1 - Sp1;

""" Discretised equation for nodal point 5 (boundary nodes) : apTp = awTw + aeTe + Su"""
aw5 = 1/dx;
ae5 = 0;
Su5 =n2*dx*Ta;
Sp5 = -n2*dx
ap5 = aw5 + ae5 - Sp5;



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
""" Exact Solution : (T-Ta)/(Tb-Ta) = cosh(n.(L-x))/cosh(nL) """
X = np.linspace(0,L,n+2);
Y = [[TbD]] + T.tolist() + [[TaD]];
x = np.linspace(0,L,50);
y = TaD + (TbD-TaD)*(np.cosh(n*(L-x))/np.cosh(n*L));


plt.grid()
plt.title('Comparaison of the Num. res. w/ the Analytical sol.')
plt.xlabel('Distance [m]')
plt.ylabel('Temperature [°C]')
plt.plot(X,Y,'r+',x,y,'b--');
plt.legend(['Numerical','Analytical']);
plt.show()

