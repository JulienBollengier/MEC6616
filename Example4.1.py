# -*- coding: utf-8 -*-
""" LAP1 """
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt



""" Example 4.1 """
""" Source-free heat conduction in an insulated rod """



""" Function to calculate the temperature and the distance where the temperature is calculated"""
""" Explanation of the different parameters :
# n is the number of points
# L is the distance of the rod en meter [m]
# k is the Thermal conductivity [W/m.K]
# A is the cross-section of the rod [m^2]
# Ta and Tb are the temperature respectively at the begining and the end of the rod"""


def Temperature(n,L,k,A,Ta,Tb):
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
    t = solve(M,S)-273.15;
    T = [[Ta-273.15]] + t.tolist() + [[Tb-273.15]];
    return T,X






""" Resolution of the exercise """
""" 
# Display the numerical solution and the analytical solution on the same graph
# Find the error between the analytical and the numerical solution
# Display the evolution ov the error with the number of points (n)
"""

""" Constants """
L = 0.5; # Rod distance in meters [m]
n = [5,10,20]; # Number of nodes [nb] """
k = 1000; # Thermal conductivity [W/m.K]
A = 10e-3; # Cross-section Area [m^2]
Ta = 100+273.15; # Temperature at the beginning in degres [K]
Tb = 500+273.15; # Temperature at the end in degres [K]


""" Numerical solution """
T = []
X = []
for i in range(len(n)):
    Tt,Xx = Temperature(n[i],L,k,A,Ta,Tb)
    T.append(Tt)
    X.append(Xx)
    


""" Anaytical solution """
""" Exact Solution : Ttheo = 800x + 100 """
x = np.linspace(0,L,50);
y = 800*x + 100;

""" Errors """
N = np.zeros(len(n))
for i in range(len(n)):
    Ttheo=[]
    S=0
    for j in range(len(X[i])):
        Ttheo.append(800*X[i][j]+100)
        S = S + abs(Ttheo[j] - T[i][j])
        N[i] = 1/n[i]*S

P = np.log(N[-1]/N[0])/np.log(n[0]/n[-1]);
    

""" Graphics Display """
plt.grid()
plt.title('Comparaison of the Num. res. w/ the Analytical sol.');
plt.xlabel('Distance [m]');
plt.ylabel('Temperature [°C]');
plt.plot(X[0],T[0],'r+',x,y,'b--');
plt.legend(['Numerical','Analytical']);
plt.show()


plt.grid()
plt.title('Comparaison of the error between Num. res. & the Analytical sol.');
plt.xlabel('ln(n)');
plt.ylabel('ln(E)');
plt.plot(np.log(np.array(n)),np.log(N));
plt.legend(['Ordre de convergence P = {}'.format(P) ]);
plt.show()