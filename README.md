**CONVECTION - DIFFUSION 1D**

On désire calculer la solution numérique de l’équation de convection-diffusion en 1D : $\frac{d}{dx}(\rho u \phi)  = \frac{d}{dx}\left(\Gamma \frac{d \phi}{dx}\right)$ 

On calculera la solution numérique à l’aide des schémas Central et Upwind, tels que présentés par Versteeg aux sections 5.3 et 5.6.
1. Trouver les coefficients $a_W$, $a_E$, $a_P$, $S_u$ et $S_p$ pour un cellule $P$ pour les cas de figures suivants :

- Cellule au centre, schéma centré et upwind
- Cellule en frontière gauche, C.L. Dirichlet, schéma centréc.
- Cellule en frontière gauche, C.L. Dirichlet, schéma upwind (u positif)
- Cellule en frontière droite, C.L. Dirichlet, schéma upwind (u positif)
- Cellule en frontière droite, C.L. Dirichlet, schéma centré 

2. Écrire un programme Python qui forme le système matriciel et qui le résout pour différentes valeurs des paramètres pour des conditions limites de Dirichlet aux deux extrémités. Reproduire les exemples 5.1 et 5.2 du livre de Versteeg et Malalasekera.
3. Déterminez l’ordre de convergence observé des méthodes Upwind et Centré

**DÉPOT SUR MOODLE**

Déposer votre programme Python sur MOODLE avant le Lundi 16 septembre 23h55. Je vais exécuter le programme et je vais vérifier que :
- Votre programme fonctionne tel qu’attendu
- Votre programme est facile à comprendre
- Le programme trace le graphe de T(x) numérique et analytique
- Le programme trace le graphe ln(E) vs ln(1/Nx)
- Le programme calcule et affiche l’ordre de convergence observé
