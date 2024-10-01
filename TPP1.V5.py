# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 16:36:34 2024

@author: tiboe
"""

import numpy as np
import time
from meshGenerator import MeshGenerator
from meshConnectivity import MeshConnectivity
from meshPlotter import MeshPlotter
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

# Définition de la classe pour la gestion des calculs de maillage
class ThermalDiffusionSolver:
    def __init__(self, length, width, k, TA, TB, source_term, lc, iterations, dz):
        self.length = length
        self.width = width
        self.k = k
        self.TA = TA
        self.TB = TB
        self.source_term = source_term
        self.lc = lc
        self.iterations = iterations
        self.dz = dz
        
        # Génération du maillage
        self.mesh_generator = MeshGenerator()
        self.mesh_parameters = {'mesh_type': 'TRI', 'lc': self.lc}
        self.mesh_obj = self.mesh_generator.rectangle([0.0, self.length, 0.0, self.width], self.mesh_parameters)

        # Conditions aux limites
        self.bcdata = (['DIRICHLET', 0], ['NEUMANN', 1], ['DIRICHLET', 2], ['NEUMANN', 3])
        
        # Connexion du maillage
        self.connectivity = MeshConnectivity(self.mesh_obj)
        self.connectivity.compute_connectivity()

        # Initialisation des matrices
        self.number_of_elements = self.mesh_obj.get_number_of_elements()
        self.A = np.zeros((self.number_of_elements, self.number_of_elements), dtype=float)
        self.B = np.zeros(self.number_of_elements, dtype=float)
        self.phi = np.zeros(self.number_of_elements, dtype=float)

    def initialize_boundary_conditions(self):
        mod = Mod(self.mesh_obj)
        phi_boundary = np.zeros(self.mesh_obj.get_number_of_boundary_faces(), dtype=float)
        phi_boundary = self.initialize_temp(phi_boundary, mod)
        return phi_boundary

    def initialize_temp(self, phi_boundary, mod):
        # Initialisation des températures aux frontières
        for condition in self.bcdata:
            if condition[0] == 'DIRICHLET':
                face_index = condition[1]
                if face_index == 0:  # Est
                    phi_boundary[face_index] = self.TA
                elif face_index == 2:  # Ouest
                    phi_boundary[face_index] = self.TB
        return phi_boundary

    def apply_boundary_conditions(self, phi_boundary):
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces()):
            self.A, self.B = self.update_boundary(self.mesh_obj, i_face, phi_boundary)

    def update_boundary(self, mesh_obj, face_index, phi_boundary):
        # Mettre à jour la matrice A et le vecteur B selon les conditions aux limites
        # Cette fonction doit être définie selon la logique que vous avez utilisée dans votre fonction originale
        return self.A, self.B

    def solve(self):
        phi_boundary = self.initialize_boundary_conditions()
        for _ in range(self.iterations):
            self.A.fill(0)
            self.B.fill(0)
            self.apply_boundary_conditions(phi_boundary)

            # Résoudre le système
            self.phi = np.linalg.solve(self.A, self.B)

    def analytical_solution(self, x):
        return ((self.TB - self.TA) / self.length + self.source_term / (2 * self.k) * (self.length - x)) * x + self.TA

    def compute_error(self, x):
        analytical_sol = self.analytical_solution(x)
        error = np.abs(self.phi - analytical_sol)
        return np.linalg.norm(error)

    def plot_results(self):
        x_values = np.linspace(0, self.length, self.number_of_elements)
        analytical_sol = self.analytical_solution(x_values)

        plt.figure()
        plt.plot(x_values, self.phi, label='Solution Numérique', marker='o')
        plt.plot(x_values, analytical_sol, label='Solution Analytique', linestyle='--')
        plt.xlabel('Position (m)')
        plt.ylabel('Température (°C)')
        plt.title('Comparaison des Solutions Numérique et Analytique')
        plt.legend()
        plt.grid()
        plt.show()

# Paramètres du problème
gamma = 0.5  # Conductivité thermique (W/m.K)
length = 0.02  # Longueur du rectangle (m)
TA = 100  # Température en A (°C)
TB = 200  # Température en B (°C)
source_term = 1000  # Terme source (kW/m³)
lc = 0.004  # Longueur caractéristique
iterations = 20
dz = 1  # Distance pour le terme source

# Exécution du solveur
solver = ThermalDiffusionSolver(length, width=length, k=gamma, TA=TA, TB=TB, source_term=source_term, lc=lc, iterations=iterations, dz=dz)
solver.solve()

# Calcul de l'erreur
x_values = np.linspace(0, length, solver.number_of_elements)
error_norm = solver.compute_error(x_values)
print(f"Erreur L2: {error_norm}")

# Analyse du temps de calcul
# Temps pour matrice dense
start_time = time.time()
phi_dense = np.linalg.solve(solver.A, solver.B)
dense_time = time.time() - start_time

# Temps pour matrice creuse
A_sparse = csr_matrix(solver.A)
start_time = time.time()
phi_sparse = spsolve(A_sparse, solver.B)
sparse_time = time.time() - start_time

print(f"Temps pour matrice dense: {dense_time:.4f} s")
print(f"Temps pour matrice creuse: {sparse_time:.4f} s")

# Visualisation des résultats
solver.plot_results()
