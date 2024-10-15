# Bînome :
n°etu : 28627431
n°etu : 28602691

# Matrice des distances entre échantillons :

|                 | GCA_000005845.2 | GCA_000008865.2 | GCA_000013265.1 | GCA_000069965.1 | GCA_030271835.1 |
|-----------------|-----------------|-----------------|-----------------|-----------------|-----------------|
| GCA_000005845.2 | 1.00000         | 0.442346001     | 0.346019362     | 0.002578191     | 0.002576698     |
| GCA_000008865.2 | 0.442346001     | 1.00000         | 0.312741620     | 0.002340611     | 0.002341023     |
| GCA_000013265.1 | 0.346019362     | 0.312741620     | 1.00000         | 0.002457800     | 0.002464379     |
| GCA_000069965.1 | 0.002578191     | 0.002340611     | 0.002457800     | 1.00000         | 0.031276364     |
| GCA_030271835.1 | 0.002576698     | 0.002341023     | 0.002464379     | 0.031276364     | 1.00000         |

# Commentaires :

Cette matrice permets de tirer des enseignements différents,
- Les valeurs proches de 1 indiquent une forte similarité entre les séquences, et inversement.
- En observant les valeurs, on peut identifier des groupes de séquences qui sont plus similaires entre elles qu'avec les autres. Ces groupes pourraient correspondre à différentes espèces, souches, ou familles de gènes.
- On peut estimer la distance évolutive entre les séquences, en comparant les distances entre différentes paires de séquences, on peut obtenir des informations sur les relations phylogénétiques entre les organismes.

# Description des méthodes :

On cherche à calculer la valeur jaccard entre plusieurs séquences. 
Pour ce faire, on construit un dictionnaire contenant les kmer des séquences qu'on compare pour compter les unions et les intersections.
Les kmer sont encodés en binaire stockés sous forme d'int. Nous comptons le reverse complement des kmer comme le même kmer, on enregistre le kmer alphabetiquement plus petit.

Les descriptions et utilisations des fonctions sont spécifiées dans les fichiers .py.