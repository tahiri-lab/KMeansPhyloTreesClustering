//
//  Super-trees.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 14/06/2022.
//

#include "Super-trees.hpp"
#include "hgt_int.hpp"
#include "K-means.hpp"

/**
 * Construit la matrice de distances Robinson-Foulds entre tous les arbres
 * du tableau mesTrees, puis lance l'algorithme K-means sur cette matrice.
 *
 * Étapes :
 *   1. Calcul du nombre de feuilles de chaque arbre via main_hgt
 *   2. Remplissage de la matrice symétrique Matrice_RF (distance RF par paire)
 *   3. Remplacement des distances négatives par : moyenne_RF + |distance|
 *      (pénalisation des paires avec peu d'espèces communes, pondérée par alpha)
 *   4. Appel de main_kmeans si le nombre d'arbres est > 3
 *   5. Libération de la mémoire et retour du temps d'exécution
 *
 * Paramètres :
 *   argv       - arguments transmis à main_kmeans
 *   tabIndices - indices de partition des arbres
 *   mesTrees   - vecteur des arbres au format Newick
 *   isBH       - true = indice Ball-Hall, false = Calinski-Harabasz
 *   alpha      - paramètre de pénalisation du chevauchement des espèces (∈ [0,1])
 *   kmin       - nombre minimal de clusters
 *   kmax       - nombre maximal de clusters
 *
 * Retourne 0 en cas de succès.
 */
int main_consense(char **argv, vector<int> tabIndices, vector <string> mesTrees, bool isBH, double alpha, int kmin, int kmax){
    time_t tbegin,tend;
    double texec=0.0;

    tbegin=time(NULL);

    double **Matrice_RF;
    double **Ww;
    double **n_identique;

    double *distances = new double[6];
    int *n_leaves = new int[mesTrees.size()+1];
    string tree;
    string tree1;
    string tree2;

    for (int j=0; j<4; j++){
        distances[j]=0.0;
    }

    /* distances[5] correspond à la valeur de alpha */
    distances[5] = alpha;

    /* Création de la matrice carrée symétrique (N×N) */
    Matrice_RF= new double*[mesTrees.size()];
    Ww= new double*[mesTrees.size()];
    n_identique= new double*[mesTrees.size()];

    for(int lineDist=0;lineDist<mesTrees.size();lineDist++){
        Matrice_RF[lineDist]= new double[mesTrees.size()];
        Ww[lineDist]= new double[mesTrees.size()];
        n_identique[lineDist]= new double[mesTrees.size()];
    }

    /* Remplissage de la matrice RF : distance RF entre chaque paire d'arbres */
    for(int line=0;line<mesTrees.size();line++){
        tree=mesTrees[line];
        main_hgt(tree,tree,distances);
        n_leaves[line] = distances[4];
    }

    for(int line=0;line<(mesTrees.size()-1);line++){
        /* Diagonale à zéro */
        Matrice_RF[line][line]=0.0;

        for(int column=(line+1);column<mesTrees.size();column++){
            tree1=mesTrees[line];
            tree2=mesTrees[column];

            /* Calcul de la distance RF entre tree1 et tree2 */
            main_hgt(tree1,tree2,distances);
            Matrice_RF[line][column]=distances[0];

            /* Nombre d'espèces communes entre tree1 et tree2 */
            n_identique[line][column]=distances[3];

            /* Symétrie de la matrice */
            n_identique[column][line]=distances[3];

            /* Remplissage symétrique sans recalcul */
            Matrice_RF[column][line]=Matrice_RF[line][column];
        }
    }

    /* Calcul de la distance RF moyenne sur les paires valides (RF >= 0) */
    double avg_RF = 0.0;
    int nb_avg_RF = 0;
    for(int line=0;line<mesTrees.size()-1;line++){
        for(int colonne=(line+1);colonne<mesTrees.size();colonne++){
            if(Matrice_RF[line][colonne]>=0){
                avg_RF += Matrice_RF[line][colonne];
                nb_avg_RF += 1;
            }
        }
    }

    avg_RF = avg_RF/(nb_avg_RF*1.0);

    /* RF = Moy(RF) + Alpha*|RF| pour les paires avec trop peu d'espèces communes */
    for(int line=0;line<mesTrees.size()-1;line++){
        for(int colonne=(line+1);colonne<mesTrees.size();colonne++){
            if(Matrice_RF[line][colonne]<0){
                Matrice_RF[line][colonne] = avg_RF + abs(Matrice_RF[line][colonne]);
                Matrice_RF[colonne][line] = Matrice_RF[line][colonne];
            }
        }
    }

    /* Initialisation de la matrice de poids Ww à 1.0 */
    for (int i=0; i<mesTrees.size(); i++){
        for (int j=0; j<mesTrees.size(); j++){
            Ww[i][j]=1.0;
        }
    }

    /* Appel de l'algorithme K-means si le nombre d'arbres est suffisant */
    if(mesTrees.size()>3){
        main_kmeans(argv,mesTrees,Matrice_RF,tabIndices,isBH,kmin,kmax);
    }

    /* Libération de la mémoire */
    for (int i=0;i<mesTrees.size();i++){
        delete [] Matrice_RF[i];
        delete [] Ww[i];
        delete [] n_identique[i];
    }
    delete [] n_leaves;
    delete [] Matrice_RF;
    delete [] Ww;
    delete [] n_identique;
    delete [] distances;

    tend=time(NULL);
    texec=difftime(tend,tbegin);
    (void)texec;

    return 0;
}
