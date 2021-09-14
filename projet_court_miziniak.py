import sys
import pandas as pd
import math
from scipy.spatial import distance_matrix


# 1)Parsing du fichier pdb - extraction des coordonnees x, y et z 
def coordonnees(file):
    """
    Fonction de parsing du fichier pdb qui extrait les coordonnées des atomes.
    """
    coord = []
    with open(file, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"): # ouverture du fichier pdb et lecture à partir des atomes
                coordonnes = {}
                if str(line[76:78].strip()) != "H": # enlève les atomes d'hydrogène
                    coordonnes["atom"] = str(line[76:78].strip())
                    coordonnes["residus"] = str(line[17:20].strip())
                    coordonnes["x"] = float(line[30:39].strip())
                    coordonnes["y"] = float(line[38:47].strip())
                    coordonnes["z"] = float(line[46:55].strip())
                    coord.append(coordonnes)
    return pd.DataFrame(coord)



#2 Fibonacci

def fibonacci_sphere(samples=1): 
    
    """
    Fonction qui permet de crée un nuage de points formant une sphere basée sur l’algorithme de Saff et
    Kuijlaars trouvé sur internet et adapté. Source : https://prograide.com/pregunta/40642/repartir-uniformement-n-points-sur-une-sphre
    """
    
    phi = math.pi * (3. - math.sqrt(5.)) #incrementation
    points = [] 
    
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2 + positiony
        radius = VdW_radius
        theta = phi * i 
        # golden angle increment 
        x = math.cos(theta) * radius + positionx
        z = math.sin(theta) * radius + positionz
            
        points.append((x, y, z)) 

             
        
    return points


#3 Fonctions de concaténation des informations
def VdW_rayon(coordonnees_atome, VdW_radius):
    VdW_atome = []
    for i in coordonnees_atome["atom"]:
        VdW_atome.append(VdW_radius[i])
    return VdW_atome

def creation_sphere(coordonnees_atome, VdW_radius):
    atom_sphere = []
    for i in coordonnees_atome["atom"]:
        atom_sphere.append(VdW_radius[i]+1.52)
    return atom_sphere 

#4 matrice de distance

def distance() :
    matrice = {}
    matrice_tot = {}
    
    for i in range(len(list_x)):
        for j in range(len(list_x)-1):
            
            
            if j > i and abs(list_x[i] - list_x[j+1])<1.52 and abs(list_y[i] - list_y[j+1])<1.52 and abs(list_z[i] - list_z[j+1])<1.52 :
            
                matrice[i] = distance_matrix([[list_x[i]],[list_y[i]],[list_z[i]]] , [[list_x[j+1]],[list_y[j+1]],[list_z[j+1]]])
        matrice_tot[i] = matrice
            
    return matrice_tot



"""
    MAIN
"""




VdW_radius = {"H" : 1.20, "C" : 1.70, "N" : 1.55, "O" : 1.52, "F" : 1.47, "P" : 1.80, "S" : 1.80,
 "Cl" : 1.75, "Cu" : 1.40}

if len(sys.argv) != 2:
    sys.exit("Rentrez un fichier pdb en argument, exemple :'python projet_court_miziniak.py 1b0q.pdb'")
    
else:
    coordonnees_atome = coordonnees(str(sys.argv[1]))
    
#5 concaténation fichier
list_atom_vdw_df = pd.DataFrame(VdW_rayon(coordonnees_atome, VdW_radius))
coord_atom_df_concat = pd.concat([list_atom_vdw_df, coordonnees_atome], axis=1)
coord_atom_df_concat.columns= ["VdW_radius",  "Atom", "res", "x", "y", "z"]
list_atom_sphere = pd.DataFrame(creation_sphere(coordonnees_atome, VdW_radius))
data = pd.concat([list_atom_sphere, coord_atom_df_concat], axis=1)
infopdb= ["Rayon_sphere", "VdW_radius", "res", "Atom", "   x  ", "   y  ", "   z  "]


#6 liste des atomes et coordonnees

list_r = coordonnees_atome['atom'].tolist()
list_x = coordonnees_atome['x'].tolist()
list_y = coordonnees_atome['y'].tolist()
list_z = coordonnees_atome['z'].tolist()

matrice = distance()

#7 création des spheres


coor_sphere = {}


for i in range(len(list_r)):
    positionx = coordonnees_atome["x"][i]
    positiony = coordonnees_atome["y"][i]
    positionz = coordonnees_atome["z"][i]
           
    if list_r[i] =="N":
        VdW_radius = 1.55
        coor_sphere[i] = fibonacci_sphere(1000)
    elif list_r[i] =="C":
        VdW_radius = 1.70
        coor_sphere[i] = fibonacci_sphere(1000)
    elif list_r[i] =="O":
        VdW_radius = 1.52
        coor_sphere[i] = fibonacci_sphere(1000)
    elif list_r[i] =="F":
        VdW_radius = 1.47
        coor_sphere[i] = fibonacci_sphere(1000)
    elif list_r[i] =="P" or list_r[i] =="S" :
        VdW_radius = 1.80
        coor_sphere[i] = fibonacci_sphere(1000)
    elif list_r[i] =="Cl":
        VdW_radius = 1.75
        coor_sphere[i] = fibonacci_sphere(1000) 
    else:
        VdW_radius = 1.40
        coor_sphere[i] = fibonacci_sphere(1000)
        

#8 calcul de de la surface    

list_vdw = data["VdW_radius"].tolist()

print("Il y a ", len(data), " atomes dans la protéine")

surface = 0
surface_atomes = {} 

for i in range(len(data)):
    surface_atomes[i] = 4 * math.pi * data["VdW_radius"][i]
    surface = surface + (4 * math.pi * data["VdW_radius"][i])

print("la surface totale de la protéines est : " , surface, " Å²")


#9 calcul de la surface absolue accessible

flag = 0

surface_innacc = 0

surface_inaccessible = 0



for i in range(len(coor_sphere)):
    flag = 0
    for j in range(len(list_x)-1):
            
        if j > i and abs(list_x[i] - list_x[j+1])<1.52 + data["VdW_radius"][i+1] + data["VdW_radius"][i] and abs(list_y[i] - list_y[j+1])<1.52 + data["VdW_radius"][i+1] + data["VdW_radius"][i] and abs(list_z[i] - list_z[j+1])<1.52 + data["VdW_radius"][i+1] + data["VdW_radius"][i]:
            flag = flag + 1
        if flag > 0:
            surface_inaccessible = ((4 * math.pi * data["VdW_radius"][i]) / 100) * flag
            

    surface_innacc = surface_innacc + surface_inaccessible

surface_acce = surface - surface_innacc

print("la surface accessible au solvant de la protéines est : " , surface_acce, " Å²")
    
#10 calcul de la surface relative accessible

surface_relative = 100 - ((surface_innacc / surface) * 100)
print("la surface relative accessible au solvant est " , surface_relative, " %")