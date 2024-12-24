# double_pendulum


Simuler le mouvement d'un pendule double est un problème de physique intéresssant, et le faire en C++ avec un rendu graphique nécessite l'utilisation de bibliothèques spécifiques pour les calculs physiques et le rendu graphique. Nous utiliserons la bibliothèque SDL (Simple DirectMedia Layer) pour le rendu graphique et les équations de Lagrange pour les calculs physiques.

Prérequis:
Installer la bibliothèque SDL2.
Comprendre les équations de Lagrange pour les pendules doubles.

Étapes:
Définir les constantes et les paramètres initiaux.
Utiliser les équations de Lagrange pour calculer les positions des pendules.
Utiliser SDL pour dessiner les pendules sur la fenêtre.

Explications:
Constantes Physiques: Définition de la constante de gravité.
Structures: Création de la structure Pendulum pour stocker les propriétés des pendules.
Update Function: Calcul des nouvelles positions des pendules en utilisant les équations de Lagrange.
SDL Initialization: Initialisation de SDL pour créer la fenêtre et le renderer.
Drawing Function: Dessin des pendules en utilisant les coordonnées calculées.
Main Loop: Boucle principale pour mettre à jour les positions et dessiner les pendules en continu.
Pour Compiler et Exécuter:
Assurez-vous d'avoir installé SDL2 sur votre système, puis compilez le programme avec:

g++ -o double_pendulum main.cpp -lSDL2 -lm
./double_pendulum

Ce programme est une base pour simuler et visualiser le mouvement d'un pendule double en temps réel. Vous pouvez l'étendre pour ajouter des fonctionnalités supplémentaires, telles que des contrôles utilisateur pour ajuster les paramètres ou enregistrer les mouvements.

