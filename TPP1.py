# -*- coding: utf-8 -*-
"""
                -- LAP3 MEC6616 Polytechnique Montreal --

    * Présentation du LAP3 :
        + Ce Laboratoire d'Apprentissage en Programmation à pour but de créer
        des maillages et d'accéder aux informations de base sur ces maillages'
        

    * Ecrit par :
        + Bollengier Julien - 2362814
        + Tiboeuf Christopher - 2362869
    
    * Date :
        + Session d'automne 2024

"""
###############################################################################
""" DEBUT DU PROGRAMME """

# IMPORT STATEMENTS
import sympy as sp 
import numpy as np
import pyvista as pv
import pyvistaqt as pvQt
from meshGenerator import MeshGenerator
from meshConnectivity import MeshConnectivity
from meshPlotter import MeshPlotter

# Paramètres du maillage
Lx = 2.0
Ly = 1.0
mesh_parameters = {'mesh_type': 'TRI', 'lc': 0.5}

# Création du maillage
mesher = MeshGenerator()
plotter = MeshPlotter()
msh_obj = mesher.rectangle([0.0, Lx, 0.0, Ly], mesh_parameters)
conec = MeshConnectivity(msh_obj)
conec.compute_connectivity()

# Paramètres physiques
phia = 1.0
phib = 0.5
GAMMA = 0,5

# Initialisation des matrices et des vecteurs
N_faces = len(msh_obj.face_to_nodes)
DA = np.zeros((N_faces, 1))
n = np.zeros((N_faces, 2))
exi = np.zeros((N_faces, 2))
enu = np.zeros((N_faces, 2))
DXI = np.zeros((N_faces, 1))
PNKSI = np.zeros((N_faces, 1))
PKSIETA = np.zeros((N_faces, 1))
D = np.zeros((N_faces, 1))
S_D_CROSS = np.zeros((N_faces, 1))

# Calcul des propriétés géométriques et physiques
for i in range(N_faces):
    arrete = msh_obj.face_to_nodes[i]
    elements = msh_obj.face_to_elements[i]
    element_gauche = elements[0]
    element_droite = elements[1]

    N1 = arrete[0]
    N2 = arrete[1]
    
    Xa = msh_obj.node_to_xcoord[N1]
    Xb = msh_obj.node_to_xcoord[N2]
    Ya = msh_obj.node_to_ycoord[N1]
    Yb = msh_obj.node_to_ycoord[N2]

    # Coordonnées des noeuds des éléments
    noeuds_gauche = msh_obj.element_to_nodes[msh_obj.element_to_nodes_start[element_gauche]:
                                             msh_obj.element_to_nodes_start[element_gauche + 1]]
    
    # Gérer les cas où il n'y a pas d'élément à droite
    if element_droite != -1:
        noeuds_droite = msh_obj.element_to_nodes[msh_obj.element_to_nodes_start[element_droite]:
                                                 msh_obj.element_to_nodes_start[element_droite + 1]]
        # Calculer les centres de gravité pour les deux éléments
        XCG_gauche = np.mean(msh_obj.node_to_xcoord[noeuds_gauche])
        YCG_gauche = np.mean(msh_obj.node_to_ycoord[noeuds_gauche])
        XCG_droite = np.mean(msh_obj.node_to_xcoord[noeuds_droite])
        YCG_droite = np.mean(msh_obj.node_to_ycoord[noeuds_droite])
    else:
        # Si pas d'élément à droite, utiliser le centre de gravité gauche seulement
        XCG_gauche = np.mean(msh_obj.node_to_xcoord[noeuds_gauche])
        YCG_gauche = np.mean(msh_obj.node_to_ycoord[noeuds_gauche])
        XCG_droite = 0.5 * (Xa + Xb)
        YCG_droite = 0.5 * (Ya + Yb)

    # Évaluation de ΔXi et autres propriétés
    DX = Xb - Xa
    DY = Yb - Ya
    DA[i, 0] = np.sqrt(DX**2 + DY**2)

    # Extraction des valeurs scalaires
    n[i, 0] = DY / DA[i, 0]  # Utilisation de DA[i, 0] pour accéder à la valeur scalaire
    n[i, 1] = -DX / DA[i, 0]  # Idem

    DXI[i, 0] = np.sqrt((XCG_droite - XCG_gauche)**2 + (YCG_droite - YCG_gauche)**2)

    # Prévention de division par zéro
    if DXI[i, 0] != 0:
        exi[i, 0] = (XCG_droite - XCG_gauche) / DXI[i, 0]
        exi[i, 1] = (YCG_droite - YCG_gauche) / DXI[i, 0]
    else:
        exi[i, :] = 0

    enu[i, 0] = DX / DA[i, 0]
    enu[i, 1] = DY / DA[i, 0]

    PNKSI[i, 0] = (DY * (XCG_droite - XCG_gauche) / (DA[i, 0] * DXI[i, 0])) - (DX * (YCG_droite - YCG_gauche) / (DA[i, 0] * DXI[i, 0]))
    PKSIETA[i, 0] = ((XCG_droite - XCG_gauche) * (Xa - Xb)) / (DXI[i, 0] * DA[i, 0]) + ((YCG_droite - YCG_gauche) * (Ya - Yb)) / (DXI[i, 0] * DA[i, 0])
    
    # Prévention de division par zéro pour D
    D[i, 0] = 1 / (PNKSI[i, 0] + 1e-10) * GAMMA * DA[i, 0] / DXI[i, 0]  # Ajout d'une petite valeur pour éviter la division par zéro

    # Gestion du flux
    epsilon = 1e-10  # Valeur seuil pour éviter la division par zéro
    S_D_CROSS[i, 0] = np.where(np.abs(PNKSI[i, 0]) > epsilon, -GAMMA * PKSIETA[i, 0] / PNKSI[i, 0] * (phib - phia), 0.0)

# Construction de la matrice des coefficients A
N_elements = len(msh_obj.element_to_nodes) // 3  # Nombre d'éléments (triangles)
A = np.zeros((N_elements, N_elements))

# Remplissage de la matrice A avec les coefficients calculés
for i in range(N_faces):
    element_gauche = msh_obj.face_to_elements[i][0]
    element_droite = msh_obj.face_to_elements[i][1]
    
    if element_gauche != -1 and element_droite != -1:
        # Assure-toi que les indices sont valides
        if element_gauche < A.shape[0] and element_droite < A.shape[1]:
            flux = S_D_CROSS[i, 0] * DA[i, 0]
            A[element_gauche, element_droite] -= flux

# Affichage de la matrice A pour vérification
print("Matrice A :")
print(A)

# Ajoute d'autres étapes selon la nécessité, comme la résolution du système, etc.


