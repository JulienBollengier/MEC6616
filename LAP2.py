"""
                -- LAP2 MEC6616 Polytechnique Montreal --

    * Présentation du LAP2 :
        + Ce Laboratoire d'Apprentissage en Programmation à pour but de caluler
        la solution numérique de l'équation de convection-diffusion en 1D. 
        + Cette solution numérique sera obtenue grâce aux schémas Central et 
        Upwind avec une condition aux limites de Dirichlet. 
        + Cette solution sera comparée avec la solution exacte.

    * Ecrit par :
        + Bollengier Julien - 2362814
        + Tiboeuf Christopher - 2362869
    
    * Date :
        + Session d'automne 2024

"""
###############################################################################
""" DEBUT DU PROGRAMME """
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve


""" DEFINITION DE LA FONCTION """
"""Cette fonction permet de calculer les coefficients aw, ae, ap, Sp et Su en 
fonction de la vitesse du fluide u, le flux de masse F et D, des propriétés de 
transport phiA et phiB et du mode : Upwind ou Central."""

def coeff(u, F, D, phiA, phiB, mode):
# SCHÉMA CENTRÉ, CELLULE AU CENTRE
    if mode == "centre":
        aw = D + F/2;
        ae = D - F/2;
        Sp = 0;
        Su = 0;
        ap = aw + ae - Sp;

# SCHÉMA UPWIND, CELLULE AU CENTRE
    if mode == "upwind":
        aw = D + F;
        ae = D;
        Sp = 0;
        Su = 0;
        ap = aw + ae - Sp;
    
# SCHÉMA CENTRÉ C.L DIRICHLET GAUCHE
    if mode == "centre_gauche":
        ae = D - F;
        aw = 0;
        Su = (2*D+F)*phiA;
        Sp = -(2*D+F)
        ap = aw + ae - Sp;

# SCHÉMA CENTRÉ C.L DIRICHLET DROITE
    if mode == "centre_droite":
        ae = 0;
        aw = D + F/2;
        Sp = -(2*D-F)
        Su = (2*D-F)*phiB
        ap = aw + ae - Sp;
        
# SCHÉMA UPWIND C.L DIRICHLET GAUCHE, u > 0
    if mode == "upwind_gauche":
        ae = D;
        aw = 0;
        Sp = -(2*D + F);
        Su = (2*D+F)*phiA;
        ap = aw + ae - Sp;
        
# SCHÉMA UPWIND C.L DIRICHLET DROITE, u > 0
    if mode == "upwind_droite":
        aw = D+F;
        ae = 0;
        Sp = -2*D;
        Su = 2*D*phiB;
        ap = aw + ae - Sp; # Si (Fe - Fa) = 0 <=> Fe = Fa.

    return aw, ae, ap, Su, Sp

""" DEFINITION DE LA MATRICE """
# Constantes :
L = 1; # Longueure du domaine en mètre 
n = 5; # Nombre de points de discrétisation
dx = L/n; # Longueure d'une section en mètre

u = 0.1; #vitesse du fluide en m/s
F = 0.1; 
D = 0.5;
phiA = 1;
phiB = 0;



i = 0;
j = dx/2;
X = np.zeros((n+2));
M = np.zeros((n,n));
S = np.zeros((n,1));


while i < n :
    if i == 0 :
        [aw, ae, ap, Su, Sp] = coeff(u, F, D, phiA, phiB, "centre_gauche");
        M[i,i] = ap
        M[i,i+1] = -ae
        S[i,:] = Su
        X[i+1] = j
    if i == n-1 :
        [aw, ae, ap, Su, Sp] = coeff(u, F, D, phiA, phiB, "centre_droite");
        M[i,i] = ap
        M[i,i-1] = -aw
        S[i,:] = Su
        X[i+2] = j + dx/2
    if i > 0 and i < n-1:
        [aw, ae, ap, Su, Sp] = coeff(u, F, D, phiA, phiB, "centre");
        M[i,i] = ap
        M[i,i+1] = -ae
        M[i,i-1] = -aw
        S[i,:] = Su
        X[i+1] = j
        X[i+2] = j + dx
    i = i+1
    j = j + dx

SOL = solve(M,S);
Y = [[phiA]] + SOL.tolist() + [[phiB]];



""" Exact Solution : Ttheo = 800x + 100 """
x = np.linspace(0,L,50);
y = (2.7183-np.exp(x))/1.7183;
#y = 1 + (1-np.exp(25*x))/7.20e10;

""" Errors """
N = np.zeros((n+2))
phitheo=[]
k=0
for j in range(len(X)):
    phitheo.append((2.7183-np.exp(X[j]))/1.7183)
    k = k + abs(phitheo[j] - Y[j])
    N[j] =1/n*k;

#P = np.log(N[-1]/N[0])/np.log(n[0]/n[-1]);





""" Graphics Display """
plt.grid()
plt.title('Comparaison of the Num. res. w/ the Analytical sol.');
plt.xlabel('Distance [m]');
plt.ylabel('Phi');
plt.plot(X,Y,'r+',x,y,'b--');
plt.legend(['Numerical','Analytical']);
plt.show()








