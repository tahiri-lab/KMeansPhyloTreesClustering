//
//  fonctions.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

/* =========================================================================
 * fonctions.cpp — Fonctions utilitaires pour la gestion des arbres et des
 * structures de données du projet K-Means Phylogenetic Trees Clustering.
 *
 * Ce fichier regroupe :
 *   - L'initialisation, l'allocation et la libération des structures InputTree
 *   - Le calcul des critères de distance entre arbres (RF, LS, BD)
 *   - La conversion de chaînes Newick en matrices de distances
 *   - La lecture des fichiers d'entrée et des paramètres
 * ========================================================================= */

#include "fonctions.hpp"
#include "utils_tree.hpp"

extern char *description;
extern int rand_bootstrap;


/**
 * Initialise une structure InputTree avec des valeurs nulles/invalides.
 *
 * Tous les pointeurs sont mis à NULL et les entiers à -1 (ou 0 pour kt),
 * ce qui permet de tester facilement si la structure a été allouée.
 */
void initInputTree(struct InputTree *aTree){
    aTree->Adjacence = NULL;
    aTree->ARETE = NULL;
    aTree->LONGUEUR = NULL;
    aTree->ADD = NULL;
    aTree->Root = -1;
    aTree->size = -1;
    aTree->SpeciesName = NULL;
    aTree->Input = NULL;
    aTree->kt = 0;
}


/**
 * Alloue la mémoire pour tous les tableaux d'une structure InputTree.
 *
 * Alloue les matrices carrées (ADD, Adjacence, Input) de taille 2n×2n,
 * la matrice de poids W de taille (n+1)×2n, le tableau ARETE, le tableau
 * LONGUEUR et le tableau de noms d'espèces SpeciesName.
 *
 * N'alloue rien si ARETE est déjà non NULL (évite les doubles allocations).
 *
 * Note : le nom "allocMemmory" (avec double 'm') est conservé tel quel
 * pour la compatibilité avec les appels existants.
 *
 * Paramètre :
 *   n - nombre d'espèces (feuilles) de l'arbre
 */
void allocMemmory(struct InputTree *aTree, int n){
    int i;

    if(aTree->ARETE == NULL){
        aTree->degre = (int*)malloc(2*n*sizeof(int));
        aTree->ADD = (double**)malloc(2*n*sizeof(double*));
        aTree->Adjacence = (double**)malloc(2*n*sizeof(double*));
        aTree->Input = (double**)malloc(2*n*sizeof(double*));
        aTree->W = (double**)malloc((n+1)*sizeof(double*));
        
        for(i=0;i<2*n;i++){
            aTree->ADD[i] = (double*)malloc(2*n*sizeof(double));
            aTree->Adjacence[i] = (double*)malloc(2*n*sizeof(double));
            aTree->Input[i] = (double*)malloc(2*n*sizeof(double));
        }

        for(i=0;i<=n;i++)
            aTree->W[i] = (double*)malloc(2*n*sizeof(double));
        
        aTree->ARETE    =(long int*)malloc(4*(2*(n))*sizeof(long int));
        aTree->LONGUEUR    =(double*)malloc((4*(n))*sizeof(double));
        aTree->SpeciesName = (char**)malloc(2*n*sizeof(char*));
        
        for(i=0;i<=n;i++)
            aTree->SpeciesName[i] = (char*)malloc(50);
    }
}


/**
 * Libère toute la mémoire allouée pour une structure InputTree complète.
 *
 * Libère les matrices ADD, Adjacence, Input (2n lignes), la matrice W
 * (n+1 lignes), les noms d'espèces (n+1 entrées), puis les tableaux
 * de premier niveau.
 *
 * Paramètre :
 *   n - nombre d'espèces utilisé lors de l'allocation
 */
void freeInputTree(struct InputTree *aTree,int n){
    int i;
    
    for(i=0;i<2*n;i++){
        free(aTree->ADD[i]);
        free(aTree->Adjacence[i]);
        free(aTree->Input[i]);
        if(i<=n){
            free(aTree->W[i]);
        }
    }
    
    for(i=0;i<=n;i++){
        free(aTree->SpeciesName[i]);
    }
    
    free(aTree->degre);
    free(aTree->ADD);
    free(aTree->Adjacence);
    free(aTree->Input);
    free(aTree->W);
    free(aTree->ARETE);
    free(aTree->LONGUEUR);
    free(aTree->SpeciesName);
}


/**
 * Libère la mémoire d'une structure InputTree réduite (sous-arbre).
 *
 * Différent de freeInputTree : la matrice Adjacence a été allouée avec
 * une marge supplémentaire (inc=10 lignes), donc sa libération porte sur
 * 2*(n+10) lignes au lieu de 2*n.
 *
 * Paramètre :
 *   n - nombre d'espèces de base (avant la marge)
 */
void freeReducedTree(struct InputTree *aTree,int n){
    int i;
    int inc = 10;   /* marge d'allocation utilisée dans CreateSubStructures */
    
    for(i=0;i<2*n;i++){
        free(aTree->ADD[i]);
        free(aTree->Input[i]);
        if(i<=n)
            free(aTree->W[i]);
    }
    
    for(i=0;i<2*(n+inc);i++){
        free(aTree->Adjacence[i]);
    }
    
    for(i=0;i<=n;i++){
        free(aTree->SpeciesName[i]);
    }
    
    free(aTree->degre);
    free(aTree->ADD);
    free(aTree->Adjacence);
    free(aTree->Input);
    free(aTree->W);
    free(aTree->ARETE);
    free(aTree->LONGUEUR);
    free(aTree->SpeciesName);
}


/**
 * Calcule les critères de distance entre deux arbres représentés par leurs
 * matrices de distances patrawise (Matrix1 et Matrix2).
 *
 * Critères calculés :
 *   - RF  : distance de Robinson-Foulds (nombre de bipartitions non partagées)
 *   - LS  : moindres carrés (somme des carrés des différences entre matrices)
 *   - BD  : distance de bipartition
 *   - QD  : distance des quartets (non calculée ici, vaut toujours 0)
 *
 * Les résultats sont stockés dans la structure aCrit.
 *
 * Paramètres :
 *   Matrix1, Matrix2 - matrices de distances des deux arbres
 *   size             - nombre d'espèces
 *   aCrit            - structure recevant les résultats
 *   L1, A1, L2, A2   - réservés pour le calcul QD, non utilisés dans cette version
 */
void computeCriteria(double ** Matrix1, double ** Matrix2, int size,struct CRITERIA *aCrit,double *L1, long int *A1,double *L2, long int *A2){

    /* L1, A1, L2, A2 non utilisés : calcul QD non implémenté */
    (void)L1; (void)A1; (void)L2; (void)A2;

    int mI,m,i,j,RF,QD = 0;
    double LS,BD=0;

    /* Distance de Robinson-Foulds */
    m=Bipartition_Table(Matrix1,aCrit->B,aCrit->PLACE,size);
    mI=Bipartition_Table(Matrix2,aCrit->BI,aCrit->PLACEI,size);
    RF = Table_Comparaison(aCrit->B,aCrit->BI,aCrit->PLACE,aCrit->PLACEI,m,mI,size);

    /* Moindres carrés : somme des carrés des différences entre matrices */
    LS = 0.0;
    for (i=1;i<=size-1;i++)
    {
        for (j=i+1;j<=size;j++){
            LS=LS + (Matrix1[i][j]-Matrix2[i][j])*(Matrix1[i][j]-Matrix2[i][j]);
            }
        }

    /* Distance de bipartition */
    BD = BipartitionDistance(aCrit->B,aCrit->BI,size);
    
    /* Distance des quartets (QD) : non calculée dans cette version.
     * Son calcul nécessite la création de deux fichiers d'arbres temporaires.
     * QD vaut -1 si non calculée, 0 ici par convention. */
    
    aCrit->LS = LS;
    aCrit->BD = BD;
    aCrit->RF = RF;
    aCrit->QD = QD;
}


/**
 * Compte le nombre d'espèces (feuilles) dans une chaîne au format Newick.
 *
 * Principe : une espèce est comptée à chaque ':' qui suit directement un nom
 * de feuille (et non une parenthèse fermante). Les longueurs de branches
 * après les nœuds internes sont ignorées grâce au témoin.
 *
 * Paramètre :
 *   newick - chaîne de caractères au format Newick (terminée par ';')
 *
 * Retourne le nombre d'espèces trouvées.
 */
int nbSpeciesNewick(string newick){
    int i=0;
    int n = 0;
    char symbol;
    char symbolOld =' ';
    int temoin =0;  /* 1 si on est après un ')' suivi d'un chiffre (longueur interne) */

    do{
        symbol = newick.at(i);
        i++;

        if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
        if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
        if(symbol==':' && temoin==1) temoin=0;
        symbolOld = symbol;
    }while(symbol != ';');

    return n;
}


/**
 * Convertit une chaîne Newick en matrice de distances patrawise (Floyd-Warshall).
 *
 * Étapes :
 *   1. Compte le nombre d'espèces (nbSpeciesNewick)
 *   2. Alloue la mémoire pour la structure (allocMemmory)
 *   3. Analyse la chaîne Newick pour extraire ARETE et LONGUEUR (lectureNewick)
 *   4. Construit la matrice d'adjacence (loadAdjacenceMatrix)
 *   5. Calcule les distances par plus court chemin (Floyd)
 *
 * La matrice résultante est stockée dans aTree->Input.
 *
 * Paramètres :
 *   newick - chaîne Newick de l'arbre
 *   aTree  - structure InputTree à remplir
 */
void newickToMatrix(string newick, struct InputTree *aTree){

    aTree->size = nbSpeciesNewick(newick);
    allocMemmory(aTree, aTree->size + 1);
    lectureNewick(newick, aTree->ARETE, aTree->LONGUEUR, aTree->SpeciesName, &aTree->kt);
    loadAdjacenceMatrix(aTree->Adjacence, aTree->ARETE, aTree->LONGUEUR, aTree->size, aTree->kt);
    Floyd(aTree->Adjacence, aTree->ADD, aTree->Input, aTree->size, aTree->kt);
}


/**
 * Lit deux arbres au format Newick, les convertit en matrices de distances
 * et les filtre pour ne garder que les espèces communes.
 *
 * Étapes :
 *   1. Convertit tree1 et tree2 en matrices via newickToMatrix
 *   2. Filtre les matrices pour ne conserver que les espèces présentes
 *      dans les deux arbres (filtrerMatrice)
 *   3. Écrit la matrice de l'arbre espèce dans tmpFile (ecrireMatrice)
 *   4. Ajoute la matrice de l'arbre gène dans tmpFile (ajouterMatriceGene)
 *
 * Paramètres :
 *   tree1, tree2       - chaînes Newick des deux arbres
 *   tmpFile            - fichier temporaire de sortie
 *   speciesTree_t      - structure InputTree pour l'arbre espèce
 *   geneTree_t         - structure InputTree pour l'arbre gène
 *   fichier_erreur     - chemin du fichier d'erreurs
 *
 * Retourne la taille finale (nombre d'espèces communes), ou -2 en cas d'erreur.
 */
int readInputFile(string tree1, string tree2, const char *tmpFile, struct InputTree *speciesTree_t, struct InputTree *geneTree_t,char *fichier_erreur){
    int finalTaille=0;

    newickToMatrix(tree1, speciesTree_t);
    newickToMatrix(tree2, geneTree_t);

    filtrerMatrice(speciesTree_t->Input, geneTree_t->Input, speciesTree_t->SpeciesName, geneTree_t->SpeciesName, speciesTree_t->size, geneTree_t->size, fichier_erreur);

    if((finalTaille=ecrireMatrice(speciesTree_t->Input, tmpFile, speciesTree_t->size, speciesTree_t->SpeciesName)) == -1)
        return -2;
    ajouterMatriceGene(geneTree_t->Input, tmpFile, geneTree_t->size, geneTree_t->SpeciesName);

    if(finalTaille < 0)
        finalTaille = 0;

    return finalTaille;
}


/**
 * Initialise et alloue une structure CRITERIA pour le calcul des distances.
 *
 * Alloue les tableaux de bipartitions B et BI (2*size-3+1 bipartitions max,
 * chacune de taille size+1), ainsi que les tableaux d'indices PLACE et PLACEI.
 * Initialise les valeurs BD, LS et RF à 0.
 *
 * Paramètres :
 *   oldCrit - structure CRITERIA à initialiser
 *   size    - nombre d'espèces
 */
void InitCriteria(struct CRITERIA * oldCrit, int size){

    int i;

    oldCrit->PLACE=(int *) malloc((2*size-3+1)*sizeof(int));
    oldCrit->PLACEI=(int *) malloc((2*size-3+1)*sizeof(int));
    oldCrit->B=(int **) malloc((2*size-3+1)*sizeof(int*));
    oldCrit->BI=(int **) malloc((2*size-3+1)*sizeof(int*));

    for (i=0;i<=2*size-3;i++)
    {
        oldCrit->B[i]=(int *) malloc((size+1)*sizeof(int));
        oldCrit->BI[i]=(int *) malloc((size+1)*sizeof(int));
    }

    oldCrit->BD=0;
    oldCrit->LS=0.0;
    oldCrit->RF=0;
}


/**
 * Libère la mémoire allouée pour une structure CRITERIA.
 *
 * Libère les tableaux de bipartitions B et BI (2*size-3 lignes),
 * puis les tableaux PLACE, PLACEI et les pointeurs B, BI.
 *
 * Paramètres :
 *   Crit - structure CRITERIA à libérer
 *   size - nombre d'espèces utilisé lors de l'allocation
 */
void FreeCriteria(struct CRITERIA * Crit,int size){

    int i;

    free(Crit->PLACE);
    free(Crit->PLACEI);
    
    for(i=0;i<=2*size-3;i++){
        free(Crit->B[i]);
        free(Crit->BI[i]);
    }
    free(Crit->B);
    free(Crit->BI);
}


/**
 * Crée les sous-structures d'un arbre (matrice d'adjacence, degrés des nœuds)
 * à partir de la matrice ADD déjà calculée.
 *
 * Étapes :
 *   1. Alloue ARETE, LONGUEUR et Adjacence si non encore alloués (avec marge de 10)
 *   2. Extrait les arêtes de l'arbre depuis ADD (Tree_edges)
 *   3. Reconstruit la matrice d'adjacence (loadAdjacenceMatrix)
 *   4. Recalcule les distances par Floyd
 *   5. Calcule le degré de chaque nœud interne
 *
 * Paramètres :
 *   aTree   - structure InputTree à compléter
 *   inc     - ignoré (forcé à 10 en interne pour la marge d'allocation)
 *   binaire - indique si l'arbre est binaire (transmis à Tree_edges)
 */
void CreateSubStructures(struct InputTree * aTree,int inc,int binaire){

    int n = aTree->size;
    int i,j;
    int kt=0;
    inc = 10;   /* marge fixe pour éviter les débordements lors du réenracinement */

    if(aTree->ARETE == NULL){
        aTree->ARETE    =(long int*)malloc(4*(2*(n+inc))*sizeof(long int));
        aTree->LONGUEUR    =(double*)malloc((4*(n+inc))*sizeof(double));
        aTree->Adjacence=(double**)malloc((2*(n+inc)+1)*sizeof(double*));
        for(i=0;i<2*(n+inc);i++)
            aTree->Adjacence[i]=(double*)malloc((2*(n+inc)+1)*sizeof(double));
    }

    kt = aTree->kt = Tree_edges(aTree->ADD,aTree->ARETE,aTree->LONGUEUR,n,binaire);

    loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,n,aTree->kt);
    Floyd(aTree->Adjacence,aTree->ADD,n,aTree->kt);

    /* Calcul du degré de chaque nœud interne */
    aTree->degre = (int*)malloc(2*(n+inc)*sizeof(int));
    for(i=1;i<=2*n-2-kt;i++){
        aTree->degre[i]=0;
        for(j=1;j<=2*n-2-kt;j++)
            if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
    }
}


/**
 * Affiche l'aide en ligne de commande du programme HGT-Detection.
 *
 * Décrit les options disponibles (critère de distance, version, enracinement,
 * scénario, nombre maximal de HGT, chemin de fichiers, etc.).
 */
void help(){
    printf("\nHGT-DETECTION version 3.0");
    printf("\nby Nadia Tahiri, Alix Boc and Vladimir Makarenkov");

    printf("\n\nUsage :\nhgt -inputfile=[inputfilename] -outputfile=[outputfilename] -criterion=[rf|ls|bd]");
    printf(" -version=[web|consol] -speciesroot=[midpoint|prompt|file] -generoot=[midpoint|prompt|file]");
    printf(" -load=[no|yes] -viewTree=[no|yes] -scenario=[unique|multiple] -nbhgt=[maxhgt] -path=[path]");

    printf("\n\ncriterion          [rf] = Robinson and Foulds distance (default)");
    printf("  \n                   [ls] = Least-Square optimization");
    printf("  \n                   [bd] = Bipartition distance");
    printf("\n\nversion            [consol] (default)");
    printf("  \n                   [web] = get result file and tree files in web ");
    printf("  \n                   format for the Trex-online web site");
    printf("\n\nspeciesroot        [midpoint] = the root is selected by the midpoint ");
    printf("  \n                   method (default)");
    printf("  \n                   [prompt] = the program ask for the root branch");
    printf("  \n                   [file] = the root is in a file called speciesroot.txt");
    printf("\n\ngeneroot           [midpoint] = the root is selected by the midpoint ");
    printf("  \n                   method (default)");
    printf("  \n                   [prompt] = the program ask for the root branch");
    printf("  \n                   [file] = the root is in a file called generoot.txt");
    printf("\n\nsubtree            [yes] = use the subtree constraint (default)");
    printf("  \n                   [no]");
    printf("\n\nscenario           [unique] (default)");
    printf("  \n                   [multiple]");
    printf("\n\nnbhgt              number max of hgts for unique and multiple scenario. ");
    printf("  \n                   Default value = 50 hgts");
    printf("\n\nload               [yes] = load the tree from files speciestree.txt ");
    printf("  \n                   genetree.txt speciesroot.txt generoot.txt");
    printf("  \n                   [no] (default)");
    printf("\n\npath               path without \"/\" at the end. Default path is \".\"");
    printf("  \n                   the path will be add to the file. ex : path/input.txt");

    printf("\n\nExample. \nCompute hgt-detection with default parameters : ./hgt -inputfile=input.txt\n\n");
}


/**
 * Vérifie si un fichier existe et est accessible en lecture.
 *
 * Retourne true si le fichier peut être ouvert, false sinon.
 *
 * Paramètre :
 *   filename - chemin du fichier à tester
 */
bool file_exists(const char * filename)
{
    FILE *file;

    if ((file=fopen(filename, "r"))==NULL)
    {
        return false;
    }
    else{
        fclose(file);
        return true;
    }
}


/**
 * Lit une matrice de distances depuis un fichier texte et la charge dans
 * une structure InputTree.
 *
 * Format du fichier :
 *   - Première ligne : taille n (nombre d'espèces)
 *   - Bloc arbre espèce : n lignes de la forme "NomEspece val1 val2 ... valn"
 *   - Bloc arbre gène   : idem
 *
 * Le paramètre Type indique quel bloc est retenu :
 *   - SPECIE : les noms et valeurs de l'arbre espèce sont copiés dans aTree
 *   - GENE   : les noms et valeurs de l'arbre gène sont copiés dans aTree
 *
 * La matrice de poids W est initialisée à 1.0 pour toutes les paires.
 *
 * Retourne 0 en cas de succès, -1 si le fichier ne peut pas être ouvert.
 *
 * Paramètres :
 *   Type - type de données à lire (SPECIE ou GENE)
 *   file - chemin du fichier d'entrée
 *   aTree - structure InputTree à remplir
 */
int readInput(int Type, const char *file,struct InputTree * aTree){

    int size,i,j;
    char name[50];
    double val;
    FILE * in;
    
    /* Ouverture du fichier */
    if((in = fopen(file,"r"))== NULL) {
        return -1;
    }

    /* Lecture de la taille de la matrice */
    fscanf(in,"%d",&size);

    /* Allocation de la mémoire */
    aTree->size = size;
    size++;   /* espace supplémentaire pour la racine */
    aTree->SpeciesName = (char **)malloc((size+1)*sizeof(char*));
    aTree->Input = (double**)malloc((2*size)*sizeof(double*));
    aTree->ADD = (double**)malloc((2*size)*sizeof(double*));
    aTree->W = (double**)malloc((size+1)*sizeof(double*));
    for(i=0;i<2*size;i++){
        aTree->ADD[i] = (double*)malloc((2*size)*sizeof(double));
        aTree->Input[i] = (double*)malloc((2*size)*sizeof(double));
        if(i<=size){
            aTree->SpeciesName[i] = (char*)malloc(SPECIES_NAME_LENGTH);
            aTree->W[i] = (double*)malloc((size+1)*sizeof(double));
        }
    }

    /* Initialisation de la matrice de poids à 1.0 */
    for(i=0;i<=size;i++)
        for(j=0;j<=size;j++)
            aTree->W[i][j] = 1.0;

    size--;

    /* Lecture du bloc arbre espèce */
    for(i = 1; i <= size; i++)
    {
        fscanf(in,"%s",name);
        if(Type == SPECIE) strcpy(aTree->SpeciesName[i],name);
        for( j = 1; j <= size; j++)
        {
            fscanf(in,"%lf",&val);
            if(Type == SPECIE) aTree->Input[i][j] = val;
        }
    }

    /* Lecture du bloc arbre gène */
    for(i = 1; i <= size; i++)
    {
        fscanf(in,"%s",name);
        if(Type == GENE) strcpy(aTree->SpeciesName[i],name);
        for( j = 1; j <= size; j++)
        {
            fscanf(in,"%lf",&val);
            if(Type == GENE) aTree->Input[i][j] = val;
        }
    }

    strcpy(aTree->SpeciesName[size+1],"Root");

    fclose(in);
    
    return 0;
}


/**
 * Initialise la structure Parameters avec les valeurs par défaut du programme.
 *
 * Tous les champs de chaînes sont remplis via snprintf pour garantir la
 * terminaison par '\0'. Les chemins de fichiers sont construits en
 * concaténant param.path avec les noms de fichiers standards.
 *
 * Valeurs par défaut notables :
 *   - criterion  : "bd" (distance de bipartition)
 *   - mode       : "multicheck"
 *   - nbhgt      : 100
 *   - bootmin    : 0
 *   - rand_bootstrap : 0
 *
 * Retourne 0.
 */
int readParameters(struct Parameters * param){

    char input[100];
    char output[100];
    char hgtResultFile[100];

    /* Valeurs par défaut des options principales */
    snprintf((*param).sort, sizeof((*param).sort), "yes");
    snprintf((*param).printWeb, sizeof((*param).printWeb), "yes");
    snprintf((*param).criterion, sizeof((*param).criterion), "bd");
    snprintf((*param).verbose, sizeof((*param).verbose), "no");
    snprintf((*param).mode, sizeof((*param).mode), "multicheck");
    snprintf((*param).viewtree, sizeof((*param).viewtree), "no");
    snprintf((*param).generoot, sizeof((*param).generoot), "midpoint");
    snprintf((*param).speciesroot, sizeof((*param).speciesroot), "midpoint");
    snprintf((*param).load, sizeof((*param).load), "no");
    snprintf((*param).version, sizeof((*param).version), "consol");
    snprintf((*param).multiple, sizeof((*param).multiple), "no");
    snprintf((*param).multigene, sizeof((*param).multigene), "no");
    (*param).path[0] = '\0';
    input[0] = '\0';
    snprintf(output, sizeof(output), "output.txt");
    snprintf(hgtResultFile, sizeof(hgtResultFile), "hgtresultfile.txt");
    snprintf((*param).scenario, sizeof((*param).scenario), "unique");
    snprintf((*param).subtree, sizeof((*param).subtree), "yes");
    snprintf((*param).bootstrap, sizeof((*param).bootstrap), "no");
    (*param).nbhgt = 100;
    (*param).bootmin = 0;
    snprintf((*param).special, sizeof((*param).special), "no");
    (*param).rand_bootstrap = 0;
    snprintf((*param).stepbystep, sizeof((*param).stepbystep), "no");

    /* Construction des chemins de fichiers à partir du chemin de base */
    snprintf((*param).inputfile, sizeof((*param).inputfile), "%s%s", (*param).path, input);
    snprintf((*param).input, sizeof((*param).input), "%sinput_.txt", (*param).path);
    snprintf((*param).outputfile, sizeof((*param).outputfile), "%s%s", (*param).path, output);
    snprintf((*param).results, sizeof((*param).results), "%sresults.txt", (*param).path);
    snprintf((*param).results_bouba, sizeof((*param).results_bouba), "%sresults2.txt", (*param).path);
    snprintf((*param).hgtResultFile, sizeof((*param).hgtResultFile), "%s%s", (*param).path, hgtResultFile);
    snprintf((*param).speciesTree, sizeof((*param).speciesTree), "%sspeciesTree.txt", (*param).path);
    snprintf((*param).geneTree, sizeof((*param).geneTree), "%sgeneTree.txt", (*param).path);
    snprintf((*param).speciesRootfile, sizeof((*param).speciesRootfile), "%sspeciesRoot.txt", (*param).path);
    snprintf((*param).speciesRootfileLeaves, sizeof((*param).speciesRootfileLeaves), "%sspeciesRootLeaves.txt", (*param).path);
    snprintf((*param).geneRootfile, sizeof((*param).geneRootfile), "%sgeneRoot.txt", (*param).path);
    snprintf((*param).errorFile, sizeof((*param).errorFile), "%serrorFile.txt", (*param).path);
    snprintf((*param).geneRootfileLeaves, sizeof((*param).geneRootfileLeaves), "%sgeneRootLeaves.txt", (*param).path);
    snprintf((*param).speciesTreeWeb, sizeof((*param).speciesTreeWeb), "%sspeciesTreeWeb.txt", (*param).path);
    snprintf((*param).geneTreeWeb, sizeof((*param).geneTreeWeb), "%sgeneTreeWeb.txt", (*param).path);
    snprintf((*param).outputWeb, sizeof((*param).outputWeb), "%soutputWeb.txt", (*param).path);
    snprintf((*param).noMoreHgtfile, sizeof((*param).noMoreHgtfile), "%snomorehgt.txt", (*param).path);
    snprintf((*param).prehgtfile, sizeof((*param).prehgtfile), "%sprehgt.txt", (*param).path);
    
    return 0;
}
