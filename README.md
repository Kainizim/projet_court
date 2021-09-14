Projet_court
========

Programme simple permettant de calculer la surface accessible d'un solvant sur une protéine

Le programme a été fait en plusieurs étapes :
1)	Parsing + Extraction des données du fichier pdb (coordonnées de chaque atome, nom, résidus) 

2)	Effectuer une recherche des atomes voisins avec une matrice de distance
 
3)	Créer un nuage de point uniformément sur la surface d’une sphère centrée sur chaque atome de la protéine. La création des sphères a été fait grâce à un code trouver sur internet à partir de l’algorithme de Saff et Kuijlaars

4)	Effectuer une boucle qui parcours tous les points d’une sphère afin de déterminer si il existe un point appartenant à une autre sphère à une distance inferieure au rayon de l’atome ciblé + le diamètre d’une molécule d’eau

5)	Ensuite il a été fait un calcul de la surface total de la protéine en effectuant un calcul de surface de chaque type d’atome en Å².

6)	Par la suite, un calcul de la surface accessible au solvant a été effectué à partir des points n’ayant pas d’atome voisin.

7)	Pour finir une comparaison a été effectuer avec un logiciel ayant le même objectif avec une protéine : 1b0q.pdb


Pour ce programme l'utilisation de l'algorithme de Saff et Kuijlaars a été 
necessaire pour la création du nuage de points autour d'une sphère.
Source : https://prograide.com/pregunta/40642/repartir-uniformement-n-points-sur-une-sphre


Pour lancer le programme il faut mettre le nom du fichier pdb en argument
exemple : $python projet_court_miziniak.py 1b0q.pdb





