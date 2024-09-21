**Structure de donnée 2D – Quadrilatères et Triangles**

Nous avons vu en classe les structures de données qui seront utilisées pour les maillages de quadrilatères et de triangles en 2D.
Vous pourrez utiliser les objets Python du module mesh6616 fourni pour créer des maillages et accéder aux informations de base sur ces maillages.
Avant de réaliser ce laboratoire, visionner les capsules vidéo et les PowerPoint sur la programmation orienté-objet en Python et la modélisation UML disponibles sur Moodle. Vous devrez réaliser vos implémentations en utilisant l’approche de programmation orienté-objet en Python et rendre un diagramme de classe UML avec votre code.
1. Tests de base sur les maillages

    a. Relation d’Euler:

La relation d’Euler pour les maillages en 2D s’écrit : $f - a + s = 1 – h$

où : 
* f : nombre de faces (nombre d’éléments quads ou triangles)
* a : nombre d’arêtes dans le maillage
* s : nombre de sommets (nombre de noeuds dans le maillage)
* h : nombre de trous (ex : 1 trou si on maille un rectangle contenant un cylindre) Écrivez un test unitaire qui pourra vérifier si le maillage généré respecte la relation d’Euler.

    b. Divergence d’un champ constant:
  
Implémentez un test pour démontrer le respect d’un champ nul pour la divergence d’un champ constant en utilisant la structure de donnée des arêtes pour construire la divergence. (voir les notes de cours « Tests de maillage et structure de donnée » sur Moodle.

2. Calcul du gradient par la méthode Least-Square

On a vu la méthode de calcul du gradient au centre des volumes de contrôle par l’approche Least- Squares. Implémentez cette méthode en utilisant principalement des boucles sur les faces du maillage.
Testez votre programme en débutant avec un champ linéaire. Dans ce cas, le calcul du gradient devrait être exact. Ensuite, initialisez des valeurs au centre des quadrilatères ou des triangles selon une fonction non-triviale. Reconstruire le gradient au centre des volumes de contrôle et comparer avec le gradient analytique en calculant une norme de l’erreur et l’ordre de convergence. Faites des essais pour différents types de conditions aux limites. On s’attend à un ordre 1 de convergence sur le module du gradient.
