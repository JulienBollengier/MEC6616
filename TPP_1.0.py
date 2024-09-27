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

mesher = MeshGenerator()
plotter = MeshPlotter()

# Paramètres du maillage
mesh_parameters = {'mesh_type': 'TRI', 'lc': 0.5}
Lx = 2.0
Ly = 1.0

# Créer un maillage
msh_obj = mesher.rectangle([0.0, Lx, 0.0, Ly], mesh_parameters)
conec = MeshConnectivity(msh_obj)
conec.compute_connectivity()


phia = 1
phib = 1
GAMMA = 1

DA = np.zeros((len(msh_obj.face_to_nodes),1))
n = np.zeros((len(msh_obj.face_to_nodes),2))
exi = np.zeros((len(msh_obj.face_to_nodes),2))
enu = np.zeros((len(msh_obj.face_to_nodes),2))
DXI = np.zeros((len(msh_obj.face_to_nodes),1))
PNKSI = np.zeros((len(msh_obj.face_to_nodes),1))
PKSIETA = np.zeros((len(msh_obj.face_to_nodes),1))
D = np.zeros((len(msh_obj.face_to_nodes),1))
S_D_CROSS = np.zeros((len(msh_obj.face_to_nodes),1))


for i in range(len(msh_obj.face_to_nodes)):
    # Arrête à l'étude
    arrete = msh_obj.face_to_nodes[i];
    
    # Elements à gauche et à droite de l'arrête étudiée
    elements = msh_obj.face_to_elements[i];
    element_gauche = elements[0];
    element_droite = elements[1];
    
    # Noeuds de l'arrete
    N1 = arrete[0];
    N2 = arrete[1];
    
    # Coordonnées des noeuds
    Xa = msh_obj.node_to_xcoord[N1]
    Xb = msh_obj.node_to_xcoord[N2]
    Ya = msh_obj.node_to_ycoord[N1]
    Yb = msh_obj.node_to_ycoord[N2]
    
    
    # Indice des noeuds des elements de gauche et de droite de l'arrête
    start_gauche = msh_obj.element_to_nodes_start[element_gauche];
    fin_gauche = msh_obj.element_to_nodes_start[element_gauche+1];
    noeuds_gauche = msh_obj.element_to_nodes[start_gauche:fin_gauche];
    
    if element_droite != -1 :
        start_droite = msh_obj.element_to_nodes_start[element_droite];
        fin_droite = msh_obj.element_to_nodes_start[element_droite+1];
        noeuds_droite = msh_obj.element_to_nodes[start_droite:fin_droite];
    
    # Coordonnées du centre des cellules à gauche et à droite de l'arrête
    XCG_gauche = 1/3*(msh_obj.node_to_xcoord[noeuds_gauche[0]]+msh_obj.node_to_xcoord[noeuds_gauche[1]]+msh_obj.node_to_xcoord[noeuds_gauche[2]]);
    YCG_gauche = 1/3*(msh_obj.node_to_ycoord[noeuds_gauche[0]]+msh_obj.node_to_ycoord[noeuds_gauche[1]]+msh_obj.node_to_ycoord[noeuds_gauche[2]]);
    
    if element_droite != -1 :
        XCG_droite = 1/3*(msh_obj.node_to_xcoord[noeuds_droite[0]]+msh_obj.node_to_xcoord[noeuds_droite[1]]+msh_obj.node_to_xcoord[noeuds_droite[2]]);
        YCG_droite = 1/3*(msh_obj.node_to_ycoord[noeuds_droite[0]]+msh_obj.node_to_ycoord[noeuds_droite[1]]+msh_obj.node_to_ycoord[noeuds_droite[2]]);
    
    # Coordonnées du centre de l'arrête
    if element_droite == -1 :
        XCG_droite = 1/2*(msh_obj.node_to_xcoord[N1] + msh_obj.node_to_xcoord[N2])
        YCG_droite = 1/2*(msh_obj.node_to_ycoord[N1] + msh_obj.node_to_ycoord[N2])
    
    
    # Evaluation de ΔXi, ΔYi et ΔAi.
    DX = Xb - Xa
    DY = Yb - Ya
    
    DA[i] = np.sqrt(DX**2 + DY**2)
    
    n[i,0] = DY/DA[i]
    n[i,1] = -DX/DA[i]
    
    DXI[i] = np.sqrt((XCG_droite - XCG_gauche)**2 + (YCG_droite - YCG_gauche)**2)
    
    exi[i,0] = (XCG_droite - XCG_gauche)/DXI[i]
    exi[i,1] = (YCG_droite - YCG_gauche)/DXI[i]
    
    enu[i,0] = DX/DA[i]
    enu[i,1] = DY/DA[i]

    PNKSI[i] = (DY*(XCG_droite-XCG_gauche)/(DA[i]*DXI[i])) - (DX*(YCG_droite-YCG_gauche)/(DA[i]*DXI[i]))

    PKSIETA[i] = ((XCG_droite-XCG_gauche)*(Xa-Xb))/(DXI[i]*DA[i]) + ((YCG_droite-YCG_gauche)*(Ya-Yb))/(DXI[i]*DA[i])
    
    D[i] = 1/PNKSI[i]*GAMMA*DA[i]/DXI[i]
    
    S_D_CROSS = -GAMMA*PKSIETA/PNKSI*(phib-phia)
    
