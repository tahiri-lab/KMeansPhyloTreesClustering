//
//  fonctions.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

#include "fonctions.hpp"
#include "utils_tree.hpp"

extern char *description;
extern int rand_bootstrap;

//===================================================================================
//============================ DEFINITION DES FONCTIONS =============================
//===================================================================================

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

//=====================================================
//=====================================================
//=====================================================

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
            /* if(i<=n)
                aTree->W[i] = (double*)malloc(2*n*sizeof(double)); */
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

void freeReducedTree(struct InputTree *aTree,int n){
    int i;
    int inc = 10;
    
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


//================================================================================
//==
//================================================================================


void computeCriteria(double ** Matrix1, double ** Matrix2, int size,struct CRITERIA *aCrit,double *L1, long int *A1,double *L2, long int *A2){

    int mI,m,i,j,RF,QD = 0;
    double LS,BD=0;

    //= robinson and foulds
    m=Bipartition_Table(Matrix1,aCrit->B,aCrit->PLACE,size);
    mI=Bipartition_Table(Matrix2,aCrit->BI,aCrit->PLACEI,size);
    RF = Table_Comparaison(aCrit->B,aCrit->BI,aCrit->PLACE,aCrit->PLACEI,m,mI,size);

    //= least-squares
    LS = 0.0;
    for (i=1;i<=size-1;i++)
    {
        for (j=i+1;j<=size;j++){
            LS=LS + (Matrix1[i][j]-Matrix2[i][j])*(Matrix1[i][j]-Matrix2[i][j]);
            if(LS > INFINI){
                
            }
        }
    }
    //= Bipartition Distance
    BD = BipartitionDistance(aCrit->B,aCrit->BI,size);
    
    //= Quartet Distance
    //= le calcul de QD necessite la creation de 2 fichiers d'arbres t1 et t2
    //= sinon la distance sera de -1 (distance non calculee)
    
    aCrit->LS = LS;
    aCrit->BD = BD;
    aCrit->RF = RF;
    aCrit->QD = QD;

}

//=================================================================
//=================================================================
//=================================================================

int nbSpeciesNewick(string newick){
    int i=0;
    int n = 0;
    char symbol;
    char symbolOld =' ';
    int temoin =0;

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


//===================================================================================
//===================================================================================
//===================================================================================


void newickToMatrix(string newick, struct InputTree *aTree){

    int pos_racine=-1;


    aTree->size = nbSpeciesNewick(newick);
    allocMemmory(aTree,aTree->size+1);
    pos_racine = lectureNewick(newick,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,&aTree->kt);

    {
        loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size,aTree->kt);
        Floyd(aTree->Adjacence,aTree->ADD,aTree->Input,aTree->size,aTree->kt); //1ere fois
    }
}

//===================================================================================
//===================================================================================
//===================================================================================

/*
 * gets two trees with a Newick r
 * **/
int readInputFile(string tree1, string tree2, const char *tmpFile, struct InputTree *speciesTree_t, struct InputTree *geneTree_t,char *fichier_erreur){
    string newick;
    int finalTaille=0;
    newick = tree1;
    newickToMatrix(newick,speciesTree_t);

    newick = tree2;
    newickToMatrix(newick,geneTree_t);

    
     filtrerMatrice(speciesTree_t->Input,geneTree_t->Input,speciesTree_t->SpeciesName,geneTree_t->SpeciesName,speciesTree_t->size,geneTree_t->size,fichier_erreur);

     if((finalTaille=ecrireMatrice(speciesTree_t->Input,tmpFile,speciesTree_t->size,speciesTree_t->SpeciesName)) == -1)
         return -2;
     ajouterMatriceGene(geneTree_t->Input,tmpFile,geneTree_t->size,geneTree_t->SpeciesName);
    
    if(finalTaille<0){
        finalTaille = 0;
    }
    
    return finalTaille;
}

//===================================================================================
//===================================================================================
//===================================================================================

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

//===================================================================================
//===================================================================================
//===================================================================================

void CreateSubStructures(struct InputTree * aTree,int inc,int binaire){

    int n = aTree->size;
    int i,j;
    int kt=0;
    inc = 10;

    //printf("n=%d",n);
    if(aTree->ARETE == NULL){
        aTree->ARETE    =(long int*)malloc(4*(2*(n+inc))*sizeof(long int));
        aTree->LONGUEUR    =(double*)malloc((4*(n+inc))*sizeof(double));
        aTree->Adjacence=(double**)malloc((2*(n+inc)+1)*sizeof(double*));
        for(i=0;i<2*(n+inc);i++)
            aTree->Adjacence[i]=(double*)malloc((2*(n+inc)+1)*sizeof(double));
    }

    //    printf("\nbinaire = %d",binaire);
    kt = aTree->kt = Tree_edges (aTree->ADD,aTree->ARETE,aTree->LONGUEUR,n,binaire);

    loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,n,aTree->kt);
    Floyd(aTree->Adjacence,aTree->ADD,n,aTree->kt); // 4eme fois
    //===creation de degre
    aTree->degre = (int*)malloc(2*(n+inc)*sizeof(int));
    for(i=1;i<=2*n-2-kt;i++){
        aTree->degre[i]=0;
        for(j=1;j<=2*n-2-kt;j++)
            if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
    }
}


//===================================================================================
//===================================================================================
//===================================================================================


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

//===================================================================================
//===================================================================================
//===================================================================================

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

//=================================================================
//== read the species or the gene tree in the input file
//=================================================================

int readInput(int Type, const char *file,struct InputTree * aTree){

    int size,i,j;
    char name[50];
    double val;
    FILE * in;
    
    //= ouverture du fichier
    if((in = fopen(file,"r"))== NULL) {
        return -1;
    }

    //= lecture de la taille des matrices
    fscanf(in,"%d",&size);
    //= allocation de la memoire
    aTree->size = size;
    size++; // more space for the root
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

    for(i=0;i<=size;i++)
        for(j=0;j<=size;j++)
            aTree->W[i][j] = 1.0;

    size--;
    //= reads species tree
    for(i = 1; i <= size; i++)
    {
        fscanf(in,"%s",name);
        if(Type == SPECIE) strcpy(aTree->SpeciesName[i],name);
        //std::cout << name << std::endl;
        for( j = 1; j <= size; j++)
        {
            fscanf(in,"%lf",&val);
            if(Type == SPECIE) aTree->Input[i][j] = val;
        }
    }

    //= read gene tree
    for(i = 1; i <= size; i++)
    {
        fscanf(in,"%s",name);
        if(Type == GENE) strcpy(aTree->SpeciesName[i],name);
        //std::cout << name << std::endl;
        for( j = 1; j <= size; j++)
        {
            fscanf(in,"%lf",&val);
            if(Type == GENE) aTree->Input[i][j] = val;
        }
    }

    strcpy(aTree->SpeciesName[size+1],"Root");

    //std::cout << "Last species : " << aTree->SpeciesName[size] << std::endl;
    
    fclose(in);
    
    return 0;
}

//=====================================================================
//=====================================================================
//=====================================================================

int readParameters(struct Parameters * param){

    char input[100];
    char output[100];
    char hgtResultFile[100];
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
    snprintf((*param).path, sizeof((*param).path), "");
    snprintf(input, sizeof(input), "");
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

