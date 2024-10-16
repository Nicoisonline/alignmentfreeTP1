# Bînome :
n°etu : 28627431
n°etu : 28602691

# Première approche : 

## Matrice des distances entre échantillons :

|                 | GCA_000005845.2 | GCA_000008865.2 | GCA_000013265.1 | GCA_000069965.1 | GCA_030271835.1 |
|-----------------|-----------------|-----------------|-----------------|-----------------|-----------------|
| GCA_000005845.2 | 1.00000         | 0.442346001     | 0.346019362     | 0.002578191     | 0.002576698     |
| GCA_000008865.2 | 0.442346001     | 1.00000         | 0.312741620     | 0.002340611     | 0.002341023     |
| GCA_000013265.1 | 0.346019362     | 0.312741620     | 1.00000         | 0.002457800     | 0.002464379     |
| GCA_000069965.1 | 0.002578191     | 0.002340611     | 0.002457800     | 1.00000         | 0.031276364     |
| GCA_030271835.1 | 0.002576698     | 0.002341023     | 0.002464379     | 0.031276364     | 1.00000         |

## Commentaires :

Cette matrice permets de tirer des enseignements différents,
- Les valeurs proches de 1 indiquent une forte similarité entre les séquences, et inversement.
- En observant les valeurs, on peut identifier des groupes de séquences qui sont plus similaires entre elles qu'avec les autres. Ces groupes pourraient correspondre à différentes espèces, souches, ou familles de gènes.
- On peut estimer la distance évolutive entre les séquences, en comparant les distances entre différentes paires de séquences, on peut obtenir des informations sur les relations phylogénétiques entre les organismes.

## Description des méthodes :

On cherche à calculer la valeur jaccard entre plusieurs séquences. 
Pour ce faire, on construit un dictionnaire contenant les kmer des séquences qu'on compare pour compter les unions et les intersections.
Les kmer sont encodés en binaire stockés sous forme d'int. Nous comptons le reverse complement des kmer comme le même kmer, on enregistre le kmer alphabetiquement plus petit.

Les descriptions et utilisations des fonctions sont spécifiées dans les fichiers .py.

# Deuxième approche :

## Matrice des distances entre échantillons :

|                 | GCA_000005845.2 | GCA_000008865.2 | GCA_000013265.1 | GCA_000069965.1 | GCA_030271835.1 |
|-----------------|-----------------|-----------------|-----------------|-----------------|-----------------|
| GCA_000005845.2 | 1.00000         | 0.392061281     | 0.303129074     | 0.004522613     | 0.004522613     |
| GCA_000008865.2 | 0.392061281     | 1.00000         | 0.283055198     | 0.003514056     | 0.003514056     |
| GCA_000013265.1 | 0.303129074     | 0.283055198     | 1.00000         | 0.004018081     | 0.004018081     |
| GCA_000069965.1 | 0.004522613     | 0.003514056     | 0.004018081     | 1.00000         | 0.046596858     |
| GCA_030271835.1 | 0.004522613     | 0.003514056     | 0.004018081     | 0.046596858     | 1.00000         |

## Commentaires :

On peut observer que les valeurs de la matrice sont différentes de celles de la première approche. Cela est dû à la méthode de calcul de la distance qui est différente.
Néanmoins, on observe que l'ordre de grandeur des valeurs est similaire, ce qui indique que les deux méthodes donnent des résultats cohérents.

## Description des méthodes :

Le but est de repartir du code de la semaine dernière et de l'optimiser. Pour ce faire, nous utilisons une méthode de hachage afin de déterminer l'union et l'intersection en comparant des listes de taille s des plus petits k-mers. L'objectif est de réduire le nombre d'éléments comparés. Pour ce faire, `stream_kmer` ne génère plus des k-mers encodés en int, mais utilise la fonction `xorshift`. La fonction `min_hash` permet de créer les listes de taille s des plus petits k-mers générés par `xorshift`, triés par ordre croissant. Enfin, dans la fonction `jaccard`, nous comparons maintenant les listes pour les fichiers A et B. Dans cette comparaison, la position du k-mer est importante car elle indique s'il y a une intersection et/ou une union. Nous calculons ensuite la valeur de Jaccard et construisons la matrice.

# Concernant les génomes Humain, souris et singe : 

La taille conséquente et le format des génomes humain, souris et singe ne nous permettent pas de les traiter de la même manière que les génomes bactériens. En effet, les génomes bactériens sont plus petits et plus simples, ce qui permet de les traiter plus rapidement. Pour ces génomes, il n'a pas été possible de calculer la matrice des distances, car cela prenait trop de temps, sans compter les difficultés liées aux données.