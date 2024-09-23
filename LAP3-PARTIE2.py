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

print('Rectangle : maillage non structuré avec des quadrilatères')
mesh_parameters = {'mesh_type': 'TRI', 'lc': 0.5}
mesh_obj = mesher.rectangle([0.0, 1.0, 0.0, 1.0], mesh_parameters)
conec = MeshConnectivity(mesh_obj)
conec.compute_connectivity()
plotter.plot_mesh(mesh_obj)


""" Fonction permettant de vérifier si le maillage générer respecte la relation d'Euler """
def Euler(coordonnees, mesh_parameters, mesh_type, rayon):
    # Maillage en fonction du type de maillage séléctionné
    if mesh_type == 'rectangle':
        msh_obj = mesher.rectangle(coordonnees,mesh_parameters)
        trous = 0;
    if mesh_type == 'back_step':
        msh_obj = mesher.back_step(coordonnees[0],coordonnees[1],coordonnees[2],coordonnees[3],mesh_parameters)
        trous = 0;
    if mesh_type == 'circle':
        msh_obj = mesher.circle(coordonnees, rayon, mesh_parameters)
        trous = 1;
    if mesh_type == 'quarter_annular':
        msh_obj = mesher.quarter_annular(coordonnees[0],coordonnees[1],mesh_parameters)
        trous = 0;
    
    # Connectivité du maillage
    conec = MeshConnectivity(msh_obj)
    conec.compute_connectivity()
    
    # Relation d'Euler 
    nb_sommets = msh_obj.node_to_xcoord
    nb_elmts = msh_obj.element_to_nodes_start
    nb_arretes = msh_obj.face_to_elements
    Euler = (len(nb_elmts)-1) - len(nb_arretes) + len(nb_sommets) + trous
    
    # Renvoie un booléen en sortie de la fonction
    return Euler == 1


""" Fonction permettant de vérifier le respect d'un champ nul pour la divergence d'un champ constant """
def Divergence(msh_obj):
    div = 0;
    divd = 0;
    divg = 0;
    for i in range(len(msh_obj.face_to_nodes)):
        arrete = msh_obj.face_to_nodes[i];

        XN1 = msh_obj.node_to_xcoord[arrete[0]];
        YN1 = msh_obj.node_to_ycoord[arrete[0]];
        XN2 = msh_obj.node_to_xcoord[arrete[-1]];
        YN2 = msh_obj.node_to_ycoord[arrete[-1]];

        N1N2 = [XN2-XN1, YN2-YN1];
        n = [-(YN2-YN1) , (XN2-XN1)];
        
        flux = np.dot(N1N2,n)
        if arrete[-1] == -1 :
            div = div + flux
        else :
            divd = divd - flux
            divg = divg + flux
            div = divd + divg
    return div == 0
    

""" Fonction permettant de déterminer le gradient au centre des volumes par la méthode des moindres-carrés """
def Gradient(msh_obj, phi0):
    ALS = np.zeros((2,2));
    NTRI = len(msh_obj.element_to_nodes_start);
    B = np.zeros((2,1));
    for i in range(len(msh_obj.face_to_nodes)):
        # Arrête à l'étude
        arrete = msh_obj.face_to_nodes[i];
        
        # Elements à gauche et à droite de l'arrête étudiée
        element_gauche = arrete[0];
        element_droite = arrete[1];
        
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
            XCG_droite = 1/2*(msh_obj.node_to_xcoord[arrete[0]] + msh_obj.node_to_xcoord[arrete[1]])
            YCG_droite = 1/2*(msh_obj.node_to_ycoord[arrete[0]] + msh_obj.node_to_ycoord[arrete[1]])
        
        # DX et DY à gauche et à droite de l'arrête
        DX_droite = XCG_gauche - XCG_droite;
        DY_droite = YCG_gauche - YCG_droite;
            
        DX_gauche = XCG_droite - XCG_gauche;
        DY_gauche = YCG_droite - YCG_gauche;
        
        DPHI_gauche = phi0;
        

        ALS[0,0] = DX_gauche * DX_gauche;
        ALS[1,0] = DX_gauche * DY_gauche;
        ALS[0,1] = DY_gauche * DX_gauche;
        ALS[1,1] = DY_gauche * DY_gauche;
        
        B[0,:] = DX_gauche * DPHI_gauche;
        B[1,:] = DY_gauche * DPHI_gauche;
        
        ATAI = np.linalg.inv(ALS);
        
        GRAD = np.dot(ATAI,B)
        
    return GRAD
        
        
       
        