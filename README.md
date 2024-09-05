**LAP1 - Laboratoire d’apprentissage programmation 1 - Semaine 1**

- Écrire un programme Python qui résout ces problèmes par la méthode des volumes finis tel que décrite par Versteeg et Malalasekera au Chapitre 4.
Faites un effort dans la programmation pour bien structurer votre code. Le nombre de points de discrétisation devra être une variable et votre programme construira automatiquement la matrice correspondant à ce nombre de points de discrétisation. Identifier clairement une section d’entrée des données, une section de définition de la matrice pour les points intérieurs, une section de définition de la matrice pour les points en frontière, une section résolution matricielle et une section de post-traitement. Ajouter suffisamment de commentaires pour rendre la lecture du code facile.

- Ces cas étant simples, il est possible d’en calculer la solution analytique. Comparer les solutions numériques T(x) obtenues par votre programme à la solution analytique pour les trois cas. Tracer sur un graphique la variation de la température T(x) versus x et superposer la solution analytique représentée par un trait plein aux solutions numériques représentées par des symboles. Observez l’évolution de la précision de la solution numérique en fonction du nombre de points de discrétisation. Pour ce faire, calculer une norme de l’erreur sur le domaine et tracer cette erreur en fonction du nombre de points. Soignez la présentation des graphiques. Calculez l’ordre de convergence observé de votre code.

× _Déposer votre programme Python sur MOODLE avant le dimanche 8 septembre 23h55. Je vais exécuter le programme et je vais vérifier que :_
  - Votre programme fonctionne tel qu’attendu
  - Votre programme est facile à comprendre
  - Le programme trace le graphe de T(x) numérique et analytique
  - Le programme trace le graphe ln(E) vs ln(1/Nx)
  - Le programme calcule et affiche l’ordre de convergence observé

    ![Capture d’écran 2024-09-05 120109](https://github.com/user-attachments/assets/72c32191-60c2-4d26-91d8-765ebd2541b9)
