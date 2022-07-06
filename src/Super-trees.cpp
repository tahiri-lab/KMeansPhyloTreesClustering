//
//  Super-trees.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 14/06/2022.
//

#include "Super-trees.hpp"
#include "hgt_int.hpp"
#include "K-means.hpp"

int main_consense(char **argv, vector<int> tabIndices, vector <string> mesTrees, int intParam, double alpha, int kmin, int kmax){
     // Variables
    time_t tbegin,tend;
    double texec=0.;

    // Start timer
    tbegin=time(NULL);                // get the current calendar time
    
/*     printf ("Consensus'tree K-means partitioning\n");
    printf("Nadia Tahiri and Vladimir Makarenkov - Departement d'informatique - Universite du Quebec a Montreal\n");
    printf ("Original code by :  Pierre Legendre - Departement de sciences biologiques - Universite de Montreal.\n");
    printf ("(c) Pierre Legendre, 1999\n");
 */
    //Varriables
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
    
    // Distances[5] correspond à la valeur de alpha
    distances[5] = alpha;
    
    //Création de la matrice carrée et symétrique (mesTrees.size()*mesTrees.size()) : Matrice_RF
    Matrice_RF= new double*[mesTrees.size()];
    Ww= new double*[mesTrees.size()];
    n_identique= new double*[mesTrees.size()];
    
    for(int lineDist=0;lineDist<mesTrees.size();lineDist++){
        Matrice_RF[lineDist]= new double[mesTrees.size()];
        Ww[lineDist]= new double[mesTrees.size()];
        n_identique[lineDist]= new double[mesTrees.size()];
    }
    
    
    // Remplissage de la matrice des distances RF en faisant appel à main_hgt
    // qui calcule la distance RF entre chaque paire d'arbre
    
    for(int line=0;line<mesTrees.size();line++){
        tree=mesTrees[line];
        main_hgt(tree,tree,distances);
        n_leaves[line] = distances[4];
    }
    
    for(int line=0;line<(mesTrees.size()-1);line++){
        //mettre des valeurs 0.0 pour la diagonale de la matrice RF
        Matrice_RF[line][line]=0.0;
        
        for(int column=(line+1);column<mesTrees.size();column++){
            // Affectation de deux arbres : tree1 et tree2
            tree1=mesTrees[line];
            tree2=mesTrees[column];
            
            // Appel des algorithmes des calcules des distances : RF
            main_hgt(tree1,tree2,distances);
            Matrice_RF[line][column]=distances[0];
            
            /*if(Matrice_RF[line][column]<=0.0){
                Matrice_RF[line][column]=0.0;
            }
            
            if(isnan(Matrice_RF[line][column])){
                Matrice_RF[line][column]=0.0;
            }*/

            //Recuperer le nombre d'espèces communes
            n_identique[line][column]=distances[3];
            
            //Et la symétrie
            n_identique[column][line]=distances[3];
            
            // pour remplir la symétrique de la matrice RF sans réaliser de calcul (car matrice carrée symétrique)
            Matrice_RF[column][line]=Matrice_RF[line][column];
            //printf("%lf ", Matrice_RF[line][column]);

        }
        //printf("\n");
    }
    
    double avg_RF = 0.0;
    int nb_avg_RF = 0;
    // RF = Mean(RF) + Alpha*Terme2, where Mean(RF) is the average of all RF between trees (where number of same leave is more than 3)
    for(int line=0;line<mesTrees.size()-1;line++){
        for(int colonne=(line+1);colonne<mesTrees.size();colonne++){
            if(Matrice_RF[line][colonne]>=0){
                avg_RF += Matrice_RF[line][colonne];
                nb_avg_RF += 1;
            }
        }
    }
    
    /*cout<<nb_avg_RF<<endl;
    cout<<avg_RF<<endl;
    avg_RF = avg_RF/(nb_avg_RF*1.0);
    cout<<avg_RF<<endl;*/
    
    avg_RF = avg_RF/(nb_avg_RF*1.0);
    // RF = Mean(RF) + Alpha*Terme2, where Mean(RF) is the average of all RF between trees (where number of same leave is more than 3)
    for(int line=0;line<mesTrees.size()-1;line++){
        for(int colonne=(line+1);colonne<mesTrees.size();colonne++){
            if(Matrice_RF[line][colonne]<0){
                Matrice_RF[line][colonne] = avg_RF + abs(Matrice_RF[line][colonne]);
                Matrice_RF[colonne][line] = Matrice_RF[line][colonne];
            }
        }
    }
    
    //printf("MATRICE RF\n");
    //for(int line=0;line<mesTrees.size();line++){
        //printf("\n");
        //for(int colonne=0;colonne<mesTrees.size();colonne++){
            //printf("%lf(%lf)", Matrice_RF[line][colonne]*(2*n_identique[line][colonne]-6),2*n_identique[line][colonne]-6);
        //}
    //}
    
    
    //creation de la matrice de distances RF : mat
    for (int i=0; i<mesTrees.size(); i++){
        for (int j=0; j<mesTrees.size(); j++){
            Ww[i][j]=1.0;
        }
    }
 
    //appel de l'algorithme de K-means:
    if(mesTrees.size()>3){
        main_kmeans(argv,mesTrees,Matrice_RF,n_identique,Ww,tabIndices,intParam,n_leaves,kmin,kmax);
    }
        
    //Liberation of memory
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
    // End timer
    tend=time(NULL);                // get the current calendar time

    // Compute execution time
    texec=difftime(tend,tbegin);    // tend-tbegin (result in second)
    /* cout<<"\nTEMPS D'EXECUTION "<<texec<<endl; */
    return 0;
}
