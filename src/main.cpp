//
//  main.cpp
//  testForTahiri
//
//  Created by Benjamin ALBERTELLI on 15/06/2022.
//

#include <iostream>
#include <boost/regex.hpp>
#include "structures.h"
#include <fstream>
#include <sstream>
#include "Super-trees.hpp"
#include "K-means.hpp"

using namespace std;

//===================================================================================
//===================================================================================
//===================================================================================

                          
int rand_bootstrap;
bool withConsensus = false;


//===================================================================================
//===================================================================================
//===================================================================================

//Prototypes of functions
int ExtraireDonneesLC(const char * chaine, char *champs, char * contenu);
void validation(int &intParam);
void validationAlpha(double &alpha);
void validationKmin(int intParam, int &kmin);
void presenterProgramme(void);
void Initialisation(int, char **);

//===================================================================================
//===================================================================================
//===================================================================================

int main(int nargs, char ** argv) {
    
    // Déclaration et initialisation de variables
    char champs[100];
    char contenu[100];

    Initialisation(nargs, argv);
    //nargs = 7;

    presenterProgramme();
    
    if(ExtraireDonneesLC(argv[1],champs,contenu)==1){

        //if program is launched with -tree argument
        if(strcmp("tree",champs) == 0){
            fstream fichier(argv[2]);
            int intParam = 0;
            double alpha = 0.0;
            int kmin=0;
            int kmax=0;
            if(nargs==3){
                intParam = 1;
                alpha = 1;
                kmin = 2;
            }else if (nargs==4){
                intParam = atoi(argv[3]);
                validation(intParam);
                alpha = 1;
                validationKmin(intParam,kmin);
            }else if (nargs==5){
                intParam = atoi(argv[3]);
                validation(intParam);
                alpha = atof(argv[4]);
                validationAlpha(alpha);
                validationKmin(intParam,kmin);
            }else if (nargs==6){
                intParam = atoi(argv[3]);
                validation(intParam);
                alpha = atof(argv[4]);
                validationAlpha(alpha);
                kmin = atoi(argv[5]);
                validationKmin(intParam,kmin);
            }else if (nargs==7){
                intParam = atoi(argv[3]);
                validation(intParam);
                alpha = atof(argv[4]);
                validationAlpha(alpha);
                kmin = atoi(argv[5]);
                validationKmin(intParam,kmin);
                kmax = atoi(argv[6]);
            }else{
                if(nargs > 7){
                    printf("\nbad input..\nusage:%s -tree nameFile [cluster_validity_index] [alpha] [kmin] [kmax]\n",argv[0]);
                    exit(1);
                }
            }

            vector <string> mesTrees;
            int ligne = 1;
            char ** cl2 = new char*[4];
            for (int i=0;i<4;i++){
                cl2[i] = new char[10];
            }

            strcpy(cl2[0], "*");
            strcpy(cl2[1], "?");
            strcpy(cl2[2], "?");
            strcpy(cl2[3], "?");
            vector <int> tabIndices;
            if( !fichier ){
                std::cout << "File "<<argv[2]<<" no exist."<<std::endl;
            }else{
                while( !fichier.eof()){
                    mesTrees.push_back("");//creation d'une ligne vide
                    getline(fichier, mesTrees.back());//lecture d'une ligne du fichier
                    ligne = int (mesTrees.size() - 1);//je recupere la taille du tableau (-1 pour la ligne 0)

                    if(mesTrees[ligne].empty())//si la ligne est vide
                        mesTrees.pop_back();//on la retire du tableau
                }
                tabIndices.push_back(int (mesTrees.size()));
                if (kmax>mesTrees.size()-1||kmax<1){
                    kmax = int (mesTrees.size()-1);
                }

                //call to main_consense, a consensus algorithm to classify super trees
                main_consense(cl2,tabIndices,mesTrees,intParam,alpha,kmin,kmax);

                //vider les vecteurs
                mesTrees.clear();
                tabIndices.clear();
            }

        //if program is launched with -matrice argument
        } else if (strcmp("matrice",champs) == 0) {
            if(nargs > 7){
                printf("\nbad input..\nusage:%s {-matrice} nameFile [cluster_validity_index] [alpha] [kmin] [kmax]\n",argv[0]);
                exit(1);
            }
            fstream fichier(argv[2]);
            int intParam = atoi(argv[3]);
            int kmin = atoi(argv[5]);
            int kmax = atoi(argv[6]);
            validation(intParam);
            vector <string> mesTrees;
            char ** cl2 = new char*[4];
            for (int i=0;i<4;i++){
                cl2[i] = new char[10];
            }

            //Varriables
            double **Matrice_RF;
            double **Ww;
            double **n_identique;
            double *distances = new double[4];
            int n = 0; //taille de la matrice RF
            string contenu = "";


            strcpy(cl2[0], "*");
            strcpy(cl2[1], "?");
            strcpy(cl2[2], "?");
            strcpy(cl2[3], "?");
            vector <int> tabIndices;
            size_t pos = 0;
            string delimiter = "\n";
            string space = " ";
            string val = "";
            int ligne = 0;
            int colonne = 0;

            if( !fichier ){
                cout << "File "<<argv[2]<<" no exist."<<endl;
            }else{

                //lecture de la premiere ligne (taille de la matrice)
                getline(fichier, contenu);//lecture d'une ligne du fichier

                val = contenu.substr(0, contenu.find(delimiter));
                istringstream(val) >> n;

                Matrice_RF= new double*[n];
                Ww= new double*[n];
                n_identique= new double*[n];

                for(int lineDist=0;lineDist<n;lineDist++){
                    Matrice_RF[lineDist]= new double[n];
                    Ww[lineDist]= new double[n];
                    n_identique[lineDist]= new double[n];
                    mesTrees.push_back("");//creation d'une ligne vide
                }


                //Initialisation des matrices : Matrice_RF, Ww et n_identique
                for (int i=0; i<n; i++)
                {
                    for (int j=0; j<n; j++)
                    {
                        Matrice_RF[i][j]=0.0;
                        Ww[i][j]=1;
                        n_identique[i][j]=n;
                    }
                }

                while( !fichier.eof()){
                    getline(fichier, contenu);//lecture d'une ligne du fichier
                    colonne = 0;
                    while ((pos = contenu.find(space))!= std::string::npos) {
                        val = contenu.substr(0, pos);
                        istringstream(val) >> Matrice_RF[ligne][colonne];

                        contenu.erase(0, pos + space.length());
                        colonne++;
                    }
                    ligne++;

                }
                tabIndices.push_back(n);

                int *n_leaves = new int[mesTrees.size()+1];
                double alpha = atof(argv[4]);
                //appel de l'algorithme de K-means:
                if(mesTrees.size()>3){
                    main_kmeans(cl2,mesTrees,Matrice_RF,n_identique,Ww,tabIndices,intParam, n_leaves,kmin,kmax, alpha);
                }

                //vider les vectors
                mesTrees.clear();
                tabIndices.clear();

                //Liberation of memory
                for (int i=0;i<n;i++){
                    delete [] Matrice_RF[i];
                    delete [] Ww[i];
                    delete [] n_identique[i];
                }
                delete [] Matrice_RF;
                delete [] Ww;
                delete [] n_identique;
                delete [] distances;
                for (int i=0;i<4;i++){
                    delete [] cl2[i];
                }
                delete [] cl2;
            }
        }

    }
    printf("END OF PROGRAM!\n\n");
    return 0;
}

//===================================================================================
//===================================================================================
//===================================================================================

int ExtraireDonneesLC(const char * chaine, char *champs, char * contenu){
    int tailleChaine;

    if(chaine[0] != '-'){
        return 0;
    }else{
        tailleChaine = (int)strlen(chaine);
        for(int i=0;i<tailleChaine;i++){
            champs[i] = chaine[i+1];
        }
    }

    return 1;
}

void validation(int &intParam){

    while(intParam<0 || intParam>2){
        cout<<"Invalid Parameter. Help (CH using by default). Help"<<endl;
        cout<<"Parameter 1: CH "<<endl;
        cout<<"Parameter 2: BH "<<endl;
        cout<<"Parameter 0: Exit "<<endl;
        intParam=1;
    }

    if(intParam==0){
        cout<<"Invalid Parameter. Help (CH using by default). Help"<<endl;
        cout<<"Parameter 1: CH (Calinski-Harabasz) "<<endl;
        cout<<"Parameter 2: BH (Ball-Hall)"<<endl;
        cout<<"Parameter 0: Exit "<<endl;
        cout<<"====END OF PROGRAM!===="<<endl;
        exit(0);
    }
}

void validationAlpha(double &alpha){
    //Validation of the alpha argument
    if(alpha<0){
        alpha=0;
    }else if(alpha>1){
        alpha=1;
    }
}

void validationKmin(int intParam, int &kmin){
    //Validation of the K-min argument
    if(intParam==1 && kmin<1){
        kmin = 2;
    }
    if(intParam==2 && kmin<1){
        kmin = 1;
    }
}

void presenterProgramme(){
    //Presentation of the Programm
   printf ("\nGenerate Tree similar\n");
    printf("Program   : KMeansSuperTreeClustering - 2021\nAuthors : Benjamin Albertelli and Nadia Tahiri - Departement d'informatique - Universite de Sherbrooke\nPresentation : This program clusters phylogenetic trees using the k-means partitioning algorithm.\nThese trees may have the same or different, but mutually overlapping, sets of leaves (the multiple supertree problem).\nPhylogenetic trees must be given in the Newick format (program input).\nA partitioning of the input trees in K clusters of trees is returned as output.\nThe optimal number of clusters can be determined either by the Calinski-Harabasz (CH) or by the Ball-Hall (BH) cluster validity index adapted for tree clustering.\nA supertree can then be inferred for each cluster of trees.\nThe Robinson and Foulds topological distance is used in the objective function of K-means.\nThe list of the program parameters is specified below.\n\n\n");
}

void Initialisation(int nargs, char ** argv){
    //Allow the user to use the program without executing the example line but choosing arguments step by step
    int choice = 0;

    std::cout<<"nargs vaut : "<<nargs<<std::endl;

    if(nargs != 7){
        if(nargs < 1){
            printf("Sorry but you have to choose a program to execute. Here are the possibilities :\n1) KMPTC\nPlease write the number of your choice :\n");
            scanf("%d", &choice);
            if (choice == 1){
                argv[0] = (char *)"./KMPTC";
            }
        }
        else{
            printf("\nSorry but apparently we can't understand what you want to execute.\nFirst, choose the type of data you want to work with :\n1) tree\n2) simulation\n3) matrix\nPlease write the number of your choice :\n");
            scanf("%d", &choice);
            if (choice == 1){
                char rep[256];
                getcwd(rep, 256);
                argv[1] = (char *)"-tree";
                printf("\nNow choose your .txt file containing the data.\n");
                printf("Please write the path directory from the src folder.\n");
                printf("Note that your actual directory path is :\n");
                printf("%s", rep);
                printf("\n");
                scanf("%s", argv[2]);
            }
            printf("\nPlease choose the cluster_validity_index used in K-means :\n1) Calinski-Harabasz\n2) Ball-Hall\n");
            scanf("%s", argv[3]);
            printf("\nα is the penalty parameter for species overlap in phylogenetic trees.\n");
            printf("It must be between 0 and 1 :\n");
            scanf("%s", argv[4]);
            printf("\nAlmost ready! You have to choose the number of cluster minimum and maximum.\n");
            if(strcmp(argv[3], "1") == 0){
                printf("You chose Calisnki-Harabasz, so Kmin has to be >= 2\n");
            }
            if(strcmp(argv[3], "2") == 0){
                printf("You chose Ball-Hall, so Kmin has to be >= 1\n");
            }
            printf("Kmin : ");
            scanf("%s", argv[5]);
            printf("Kmax : ");
            scanf("%s", argv[6]);
        }
    }
}
