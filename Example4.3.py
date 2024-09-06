# -*- coding: utf-8 -*-
""" LAP1 """

""" Example 4.1 """
""" Cooling of a circular fin by means of convective heat transfer along its length. """
import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt


""" Function to calculate the temperature and the distance where the temperature is calculated"""
""" Explanation of the different parameters :
# n is the number of points
# L is the distance of the rod en meter [m]
# k is the Thermal conductivity [W/m.K]
# A is the cross-section of the rod [m^2]
# Ta and Tb are the temperature respectively at the begining and the end of the rod"""

def Temperature(n,L,n2,Ta,Tb):    
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
    T = [[TbD]] + t.tolist() + [[TaD]];
    return T,X

""" Resolution of the exercise """
""" 
# Display the numerical solution and the analytical solution on the same graph
# Find the error between the analytical and the numerical solution
# Display the evolution ov the error with the number of points (n)
"""
""" Constants """
L = 1; # Rod distance in meters [m]
n = [5,10,15]; # Number of nodes [nb] 
q = 1000e3; # Heat Generation [W/m^3]
n2 = np.array(n)**2; # n^2 = h.P/k.A [m^(-2)] 
Ta = 293.15; # Ambiant Temperature in degres [K]
Tb = 373.15; # Temperature at the end in degres [K]
TaD = Ta - 273.15; # Temperature at the beginning in degres [°C] 
TbD = Tb - 273.15; # Temperature at the end in degres [°C] 

""" Numerical solution """
T = []
X = []
for i in range(len(n)):
    Tt,Xx = Temperature(n[i],L,n2[i],Ta,Tb)
    T.append(Tt)
    X.append(Xx)

""" Anaytical solution """
""" Exact Solution : (T-Ta)/(Tb-Ta) = cosh(n.(L-x))/cosh(nL) """
x = np.linspace(0,L,50);
y = TaD + (TbD-TaD)*(np.cosh(n[0]*(L-x))/np.cosh(n[0]*L));

""" Errors """
N = np.zeros(len(n))
for i in range(len(n)):
    Ttheo=[]
    S=0
    for j in range(len(X[i])):
        Ttheo.append((TaD + (TbD-TaD)*(np.cosh(n[i]*(L-X[i][j]))/np.cosh(n[i]*L))))
        S = S + abs(Ttheo[j] - T[i][j])
        N[i] = 1/n[i]*S

P = np.log(N[-1]/N[0])/np.log(n[0]/n[-1]);

""" Graphs """
plt.grid()
plt.title('Comparaison of the Num. res. w/ the Analytical sol.')
plt.xlabel('Distance [m]')
plt.ylabel('Temperature [°C]')
plt.plot(X[0],T[0],'r+',x,y,'b--');
plt.legend(['Numerical','Analytical']);
plt.show()

plt.grid()
plt.title('Comparaison of the error between Num. res. & the Analytical sol.');
plt.xlabel('ln(1/n)');
plt.ylabel('ln(E)');
plt.plot(np.log(1/np.array(n)),np.log(N));
plt.legend(['Ordre de convergence P = {}'.format(P) ]);
plt.show()