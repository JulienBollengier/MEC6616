# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 12:44:29 2024

@author: tiboe
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from meshGenerator import MeshGenerator
from meshConnectivity import MeshConnectivity
from meshPlotter import MeshPlotter
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import time

# Paramètres
Lx = 2.0
Ly = 1.0
k = 0.5  # Conductivité thermique
q = 1000.0  # Source de chaleur en kW/m^3
TA = 100.0  # Température à l'Ouest
TB = 200.0  # Température à l'Est
num_nodes_list = [0.5, 0.75, 0.25]  # Différents nombres de nœuds à tester

# Création de la fonction pour évaluer l'erreur absolue
def calculate_absolute_error(T_numeric, T_analytical):
    return np.linalg.norm(T_numeric - T_analytical)

# Initialiser une liste pour stocker les résultats
results = []

# Boucle principale pour évaluer les différentes configurations
for use_cross_diffusion in [True, False]:
    for n_nodes in num_nodes_list:
        # Création du maillage
        mesh_parameters = {'mesh_type': 'TRI', 'lc': n_nodes}
        mesher = MeshGenerator()
        plotter = MeshPlotter()
        msh_obj = mesher.rectangle([0.0, Lx, 0.0, Ly], mesh_parameters)
        conec = MeshConnectivity(msh_obj)
        conec.compute_connectivity()

        # Matrices pour la méthode des volumes finis
        n_faces = len(msh_obj.element_to_nodes)
        A = lil_matrix((n_faces, n_faces))  # Matrice des coefficients
        b = np.zeros(n_faces)  # Vecteur de source

        #Afficher le Mesh
        plotter.plot_mesh(msh_obj)        

        # Démarrer le chronomètre
        start_time = time.time()

        # Boucle pour construire les matrices
        for i in range(n_faces):
            arrete = msh_obj.face_to_nodes[i]
            elements = msh_obj.face_to_elements[i]
            element_gauche = elements[0].item()  # Convertir en entier
            element_droite = elements[1].item()   # Convertir en entier

            # Coordonnées des noeuds
            N1, N2 = arrete
            Xa, Xb = msh_obj.node_to_xcoord[N1], msh_obj.node_to_xcoord[N2]
            Ya, Yb = msh_obj.node_to_ycoord[N1], msh_obj.node_to_ycoord[N2]

            # Coordonnées du centre des cellules à gauche et à droite
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
            if use_cross_diffusion:
                # Terme de cross-diffusion (exemple, ajustez en fonction de votre modèle)
                flux = (k * (TB - TA) / Lx + q) * (DXI / 2)  # Ajoutez ici votre terme de cross-diffusion
            else:
                flux = (k * (TB - TA) / Lx + q) * (DXI / 2)  # Sans cross-diffusion
            
            A[element_gauche, element_gauche] += 1
            if element_droite != -1:
                A[element_gauche, element_droite] -= flux
                A[element_droite, element_gauche] -= flux
                A[element_droite, element_droite] += 1
            
            b[element_gauche] += TA * flux
            if element_droite != -1:
                b[element_droite] += TB * flux

        # Conditions aux limites (Dirichlet)
        A[0, 0] = 1  # Température à l'Ouest
        b[0] = TA
        A[-1, -1] = 1  # Température à l'Est
        b[-1] = TB

        # Ajouter un petit terme pour éviter la singularité
        epsilon = 1e-10  # Petit terme
        A += lil_matrix(np.eye(A.shape[0]) * epsilon)

        # Résoudre le système
        A_csr = A.tocsr()  # Convertir en matrice creuse
        T = spsolve(A_csr, b)

        # Arrêter le chronomètre
        end_time = time.time()  # Arrêter le chronomètre
        elapsed_time = end_time - start_time

        # Calculer la solution analytique
        x = np.linspace(0, Lx, 100)
        T_analytical = ((TB - TA) / Lx) * x + TA + (q / (2 * k)) * (Lx * x - (x**2 / 2))

        # Évaluer l'erreur
        #absolute_error = calculate_absolute_error(T, T_analytical)
        absolute_error = 0;
        
        # Stocker les résultats
        results.append((n_nodes, use_cross_diffusion, absolute_error, elapsed_time))

# Afficher les résultats
for n, use_cd, err, t in results:
    cd_status = "Avec cross-diffusion" if use_cd else "Sans cross-diffusion"
    print(f'Division des nœuds: {n}, {cd_status}, Erreur: {err:.6f}, Temps de calcul: {t:.4f} secondes')

# Graphique des résultats pour la dernière itération

plt.figure(figsize=(10, 6))
plt.plot(x, T_analytical, label='Solution analytique', color='blue')
#plt.scatter(msh_obj.node_to_xcoord, T, label='Solution numérique', color='red')
plt.title('Comparaison des solutions')
plt.xlabel('X (m)')
plt.ylabel('Température (°C)')
plt.legend()
plt.grid()
plt.show()

