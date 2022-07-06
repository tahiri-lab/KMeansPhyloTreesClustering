//
//  structures.h
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

#ifndef structures_h
#define structures_h


//================ CONSTANTES ========================

#define SPECIE 1
#define GENE 2
#define SPECIES_NAME_LENGTH 50
#define epsilon_ba 0.00001
#define INFINI 999999.99
#define MaxRF 0
#define MAX_HGT 100
#define FAIL -1
#define TRUE 1
#define FALSE 0
#define TRIVIAL 4


//================= STRUCTURES =======================

struct Parameters{
    char special[10];
    char printWeb[10];
    char inputfile[200];
    char results[200];
    char results_bouba[200];
    char outputfile[200];
    char hgtResultFile[200];
    char criterion[10];
    char version[20];
    char generoot[20];
    char speciesroot[10];
    char load[4];
    char viewtree[4];
    char multiple[4];
    char path[50];
    char input[200];
    char speciesTree[200];
    char geneTree[200];
    char errorFile[200];
    char noMoreHgtfile[200];
    char speciesRootfile[200];
    char speciesRootfileLeaves[200];
    char geneRootfile[200];
    char geneRootfileLeaves[200];
    char speciesTreeWeb[200];
    char geneTreeWeb[200];
    char outputWeb[200];
    char prehgtfile[200];
    char subtree[4];
    char scenario[100];
    int nbhgt;
    char bootstrap[4];
    char multigene[4];
    char mode[50];
    char verbose[50];
    char sort[10];
    int bootmin;
    int rand_bootstrap;
    char stepbystep[10];
};

struct InputTree{
    int size;
    double ** ADD;
    double ** Input;
    char ** SpeciesName;
    int Root;
    double ** Adjacence;
    long int * ARETE;
    double * LONGUEUR;
    double **W;
    int kt;
    int *degre;
};

struct CRITERIA{

    double LS;
    double rLS;
    double BD;
    double rBD;
    int rRF;
    double diff_bd;
    int RF;
    int QD;
    int m;
    int nbHgtFound;
    int ** B;
    int ** BI;
    int * PLACE;
    int * PLACEI;
};

struct HGT{
    int valide;
    int source;
    int destination;
    int * listSource;        //== 0 give the number of elements in the array
    int * listDestination;  //== 0 give the number of elements in the array
    struct CRITERIA crit;
    int source_A;
    int source_B;
    int dest_A;
    int dest_B;
    int trivial;
    int sequence;
    int source_Ar;
    int source_Br;
    int dest_Ar;
    int dest_Br;
};

struct DescTree
{
    int nbSommet;
    double ** Matrice;
    int *Tableau;
};
                                                                         
struct TNoeud
{
    int NoNoeud;
    struct TNoeud **fils;
    int nbfils;
};

struct ReduceTrace
{
    int size;
    int *species;
    int *gene;
    int *map;
    int *speciesToGene;
};

struct TreeHGT
{
    struct HGT * tabHGT;
    struct TreeHGT * fils;
    struct TreeHGT * parent;
    int taille;
};

struct TNode{
    int NoNoeud;
    char *seq;
    TNode * gauche;
    TNode * droit;
    int nbFils;
};

#endif /* structures_h */
