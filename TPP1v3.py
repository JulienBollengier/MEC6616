# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:29:21 2024

@author: tiboe
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.interpolate import griddata
from meshGenerator import MeshGenerator
from meshConnectivity import MeshConnectivity
from meshPlotter import MeshPlotter

# Paramètres
Lx = 2.0
Ly = 2.0
k = 0.5  # Conductivité thermique
q = 1000.0  # Source de chaleur en kW/m^3
TA = 100.0  # Température à l'Ouest
TB = 200.0  # Température à l'Est
n_nodes = 4  # Nombre de nœuds sur Lx
mesh_parameters = {'mesh_type': 'TRI', 'lc': Lx/n_nodes}  # Utiliser TRI pour un maillage de triangles

# Création du maillage
mesher = MeshGenerator()
plotter = MeshPlotter()
msh_obj = mesher.rectangle([0.0, Lx, 0.0, Ly], mesh_parameters)
conec = MeshConnectivity(msh_obj)
conec.compute_connectivity()

# Matrices pour la méthode des volumes finis
n_faces = len(msh_obj.face_to_nodes)
n_elements = len(msh_obj.element_to_nodes)
A = lil_matrix((n_elements, n_elements))  # Matrice des coefficients
b = np.zeros(n_elements)  # Vecteur de source

# Démarrer le chronomètre
start_time = time.time()

# Boucle pour construire les matrices
for i in range(n_faces):
    arrete = msh_obj.face_to_nodes[i]
    elements = msh_obj.face_to_elements[i]
    element_gauche = elements[0].item()  # Convertir en entier
    element_droite = elements[1].item() if len(elements) > 1 else -1   # Convertir en entier, -1 si pas d'élément à droite

    # Coordonnées des noeuds
    N1, N2 = arrete
    Xa, Xb = msh_obj.node_to_xcoord[N1], msh_obj.node_to_xcoord[N2]
    Ya, Yb = msh_obj.node_to_ycoord[N1], msh_obj.node_to_ycoord[N2]

    # Coordonnées du centre des cellules à gauche
    XCG_gauche = np.mean(msh_obj.node_to_xcoord[msh_obj.element_to_nodes[msh_obj.element_to_nodes_start[element_gauche]:
                                                                     msh_obj.element_to_nodes_start[element_gauche + 1]]])
    YCG_gauche = np.mean(msh_obj.node_to_ycoord[msh_obj.element_to_nodes[msh_obj.element_to_nodes_start[element_gauche]:
                                                                     msh_obj.element_to_nodes_start[element_gauche + 1]]])

    if element_droite != -1:
        XCG_droite = np.mean(msh_obj.node_to_xcoord[msh_obj.element_to_nodes[msh_obj.element_to_nodes_start[element_droite]:
                                                                         msh_obj.element_to_nodes_start[element_droite + 1]]])
        YCG_droite = np.mean(msh_obj.node_to_ycoord[msh_obj.element_to_nodes[msh_obj.element_to_nodes_start[element_droite]:
                                                                         msh_obj.element_to_nodes_start[element_droite + 1]]])
    else:
        # Pour les frontières de domaine
        XCG_droite = (Xa + Xb) / 2
        YCG_droite = (Ya + Yb) / 2

    # Évaluation des longueurs et des normales
    DX = Xb - Xa
    DY = Yb - Ya
    DA = np.sqrt(DX**2 + DY**2)

    # Normalisation des vecteurs
    n = np.array([DY / DA, -DX / DA]).T
    DXI = np.sqrt((XCG_droite - XCG_gauche)**2 + (YCG_droite - YCG_gauche)**2)
    exi = (XCG_droite - XCG_gauche) / DXI if DXI != 0 else 0
    enu = (YCG_droite - YCG_gauche) / DXI if DXI != 0 else 0

    # Assemblage de la matrice A
    flux = (k * (TB - TA) / Lx + q) * (DXI / 2)  # Flux calculé
    A[element_gauche, element_gauche] += 1
    if element_droite != -1:
        A[element_gauche, element_droite] -= flux
        A[element_droite, element_gauche] -= flux
        A[element_droite, element_droite] += 1

    b[element_gauche] += TA * flux
    if element_droite != -1:
        b[element_droite] += TB * flux

A = np.zeros((n_elements, n_elements), dtype=float)
B = np.zeros(n_elements, dtype=float)
source = np.zeros(n_elements)
    
# Boundary faces
for i_face_front in range(msh_obj.get_number_of_boundary_faces()):
    A, B = fonctions.boundary(msh_obj, mod, i_face_front, bcdata, gamma, TA, TB, L, A, B, phi_frontiere, q, dz)
# Internal faces
for i_face_interne in range(msh_obj.get_number_of_boundary_faces(), mesh_obj.get_number_of_faces()):
    A, B = fonctions.interne(mesh_obj, mod, i_face_interne, gamma, A, B, GRAD, q, dz)
        
# Ajout du terme source
for i_element in range(n_elements):
        start = msh_obj.element_to_nodes_start[i_element]
        end = msh_obj.element_to_nodes_start[i_element+1]
        nodes = msh_obj.element_to_nodes[start:end]
        # print(nodes)
        source[i_element] = q*np.abs(mod.aire_element(nodes))
        # print(nodes)
        # print("aire",np.abs(mod.aire_element(nodes)))
        # print(q*np.abs(mod.aire_element(nodes)))
        B[i_element] += source[i_element]
        
phi = np.linalg.solve(A, B)
# print("\nmax diff phi-phi2",np.max(np.abs(phi2-phi)))
# phi = phi2
GRAD = ls.least_square(msh_obj, bcdata, phi, phi_frontiere, q, dz)


# Résoudre le système linéaire
try:
    phi = np.linalg.solve(A, B)
except np.linalg.LinAlgError:
    print("Matrice A singulière. Vérification des valeurs de A et B :")
    print("Matrice A :")
    print(A)
    print("Vecteur B :")
    print(B)


# Arrêter le chronomètre
end_time = time.time()  # Arrêter le chronomètre
elapsed_time = end_time - start_time
print(f'Temps de calcul : {elapsed_time:.4f} secondes')

# Calculer les coordonnées des centres des éléments
element_centers_x = []
element_centers_y = []

for element in range(n_elements):
    start_idx = msh_obj.element_to_nodes_start[element]
    end_idx = msh_obj.element_to_nodes_start[element+1] 
    nodes = msh_obj.element_to_nodes[start_idx:end_idx]
    
    center_x = np.mean(msh_obj.node_to_xcoord[nodes])
    center_y = np.mean(msh_obj.node_to_ycoord[nodes])
    
    element_centers_x.append(center_x)
    element_centers_y.append(center_y)

element_centers_x = np.array(element_centers_x)
element_centers_y = np.array(element_centers_y)

# Interpolation pour obtenir T sur les nœuds
nodes_x = msh_obj.node_to_xcoord
nodes_y = msh_obj.node_to_ycoord

# Interpolation des valeurs de température
if len(T) == len(element_centers_x) and len(T) == len(element_centers_y):
    T_interpolated = griddata((element_centers_x, element_centers_y), T, (nodes_x, nodes_y), method='linear')
else:
    print("Erreur : Les tailles ne correspondent pas pour l'interpolation !")

# Afficher les résultats
plt.figure(figsize=(10, 6))
plt.scatter(nodes_x, T_interpolated, label='Solution numérique (interpolée)', color='red')
plt.title('Distribution de température')
plt.xlabel('X (m)')
plt.ylabel('Température (°C)')
plt.legend()
plt.grid()
plt.show()
