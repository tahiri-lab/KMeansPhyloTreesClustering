//
//  K-means.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 14/06/2022.
//
// =============================================================================================================
// Program   : KMeansSuperTreeClustering - 2018
// Authors   : Nadia Tahiri and Vladimir Makarenkov (Université du Québec a Montréal)
// This program clusters phylogenetic trees using the k-means partitioning algorithm.
// These trees may have the same or different, but mutually overlapping, sets of leaves (the multiple supertree problem).
// Phylogenetic trees must be given in the Newick format (program input).
// A partitioning of the input trees in K clusters of trees is returned as output.
// The optimal number of clusters can be determined either by the Calinski-Harabasz (CH) or by the Ball-Hall (BH) cluster
// validity index adapted for tree clustering.
// A supertree can then be inferred for each cluster of trees.
// The Robinson and Foulds topological distance is used in the objective function of K-means.
// The list of the program parameters is specified in the file README.md
// =============================================================================================================

#include "K-means.hpp"
#include "hgt_int.hpp"

// =============================================================================================================
// CONSTANTS
// =============================================================================================================
const double MIN_CH_VALUE = -1e10;              // Minimum Calinski-Harabasz index value
const double MAX_W_VALUE = 1e9;                 // Maximum W clustering validity index value
const double MAX_FO_VALUE = 1e10;               // Maximum objective function value
const double INITIAL_SSE_REF = 1.0e20;          // Initial sum of squared errors reference value
const double INITIAL_MAX_CH = -1e10;            // Initial maximum CH value
const int MAX_ITERATIONS = 100;                 // Maximum K-means iterations per cluster size
const int RAND_MAX_VALUE = 32767;               // Maximum random number from rand()
const int MAX_FILENAME_LENGTH = 300;            // Maximum length for filename strings
const int MAX_PATH_LENGTH = 255;                // Maximum length for file paths
const int ROUNDING_PRECISION = 3;               // Decimal places for rounding distance values
const int CONVERGENCE_THRESHOLD_DIVISOR = 1000; // Divisor for convergence threshold check
// =============================================================================================================

FILE *Output4;


//  Parameters of the program:
//  n    = number of observations (ex. number of trees)
//  nmax = maximum number of observations (default: 10000)
//  p    = number of variables     (ex. number of variables describing each trees)
//  pmax = maximum number of variables  (default: 10000)
//  k    = number of groups (centroids)
//  k_capacity = maximum number of groups
//  currentK   = nombre courant de clusters (groupes) en cours d’évaluation
//  MAX_ITERATIONS = maximum iteration for convergeance of centroid (fixed=100)
//  Parameter (nmax=100000,pmax=10,k_capacity=10)
//  critera = (0,1,2)
//            0: C-H
//            1: logSS
//            2: Silhouette
//              3: W
//   mat       =  double[nmax+1][pmax+1] (data matrix))
//   weight         =  double vector[p] of weigth vector for variable p
//   xbar      =  vector [k][i] = mean distance of i to centroid k
//   list      =  vector[i]=k contains group assignments for the objects.
//   howmany[k]=  contains the number of objects in each group centroid (k).
//   nobest
//   var       = double [k_capacity+1]
//   vect      = double [pmax+1]
//   mean[p]      = vector of mean value not weigthed for the variable [p] double [pmax+1]
//   Dvec      = double [k_capacity+1]
//   CHr       = vector[k] of the best CH this k
//   BHr       = vector[k] of the best BH this k
//   Silr      = vector[k] of the best Silr this k
//   LogSSr    = vector[k] of the best LogSS this k
//   Wr        = vector[k] of the best W this k
//   SSEr      = double [k_capacity+1]  sum of squared error statistic (SSE) = within-group sum of squares
//   listr     = int[k_capacity+1][nmax+1];
// sx = new double*[k_capacity+1];
//  sx2 = new double*[k_capacity+1];

//   no = new int [nmax+1];
//    iordre = new int [nmax+1];
//   nobest = new int [k_capacity+1];
//   nnitr = new int [k_capacity+1];
//
// Output file for summary results
//
// Output2   = file.statistics.txt
// Output3   = file.sum.txt
// Output4   = stat.txt
// Output7    = file.groups.txt
//
// Note, we expect the input file to be a matrix in the form (below))
// were 5 is the n, and 4 is the p and the following line are the value (double))
//5    4
//0.0    0.0    0.0    1.0    0.0
//0.0    1.0    0.0    0.0    1.0
//1.0    1.0    1.0    0.0    1.0
//0.0    0.0    0.0    0.0    1.0
//
// IMPORTANT: See also the variables below

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

struct context {
    int treeAmount;
    int k1 = 0;
    int k2 = 0;
    int nmax;
    int pmax;
    int k_capacity;
    vector<int> Strouve;
    vector<vector<int>> listr;
};

/*
Cette structure fais l'itilialisation de nombreuses variables qui seront utilisées dans main_kmeans.
Quand la transformation de main_kmeans en plein d'appel de fonctions sera terminée,
on aura pas besoin de précéder la plupart des variables avec un le nom de la structure car ceci cera fait dans la signature de la fonction.
Pour l'instant, je dois faire des truc stupide comme "treeAmount = ctx.treeAmount" dans main_kmeans pour vérifier que mon refactoring ne change
pas le comportement du code sans avoir à replacer tout les "treeAmount" par "ctx.treeAmount".
*/
context initialisation(vector <string> monTableau) {
    context ctx;
    ctx.treeAmount = int(monTableau.size()); //quantity of initial tree
    ctx.nmax = ctx.treeAmount;
    ctx.pmax = ctx.treeAmount;
    ctx.k_capacity = ctx.treeAmount;

    ctx.Strouve.reserve(ctx.treeAmount);
    for (int i=0;i<ctx.treeAmount;i++) {
        ctx.Strouve[i] = 0;
    }

    vector<int> temp1_listr((ctx.nmax+1),1);
    vector<vector<int>> temp2_listr((ctx.k_capacity+1),(temp1_listr));
    ctx.listr = temp2_listr;

    return ctx;
}

int main_kmeans(char **argv, vector <string> monTableau, double ** mat, vector<int> tabIndices, bool isBH, int k_min, int k_max){
    //*****************Define variables******************************************//
    // Variables
    time_t tbegin2 = time(NULL);

    context ctx = initialisation(monTableau);

    double CHr_max = INITIAL_MAX_CH;
    int CHr_group = 0;

    double W_min = MAX_W_VALUE;
    double W_max = MIN_CH_VALUE;
    int W_group = 0;
    double FO_new = MAX_FO_VALUE;

    const int treeAmount = ctx.treeAmount;
    int currentK=0; //nombre courant de clusters (groupes) en cours d’évaluation
    bool debug=false;
    int k1=ctx.k1, k2=ctx.k2;

    int random_number=100; //--Fixed random number
    int iassign=2;  // 1 equal, 2 random
    int nran=100;  //--Number of Random start VM

    int nmax=ctx.nmax;    //--Maximum number of object -Parameter (nmax=10000,pmax=250,k_capacity=100)
    int pmax=ctx.pmax;      //--Maximum data point (variable))
    int k_capacity=ctx.k_capacity;      // Maximum number of groups

    //char *criteria = argv[0];
    // ------------------------------------------------------------
    // Correction de la securite memoire
    // On ne doit pas modifier argv[0] car il pointe vers une zone mémoire qui peut être en lecture seule.
    // On crée donc un buffer local modifiable.
    // ------------------------------------------------------------
    char criteria_buf[8] = {0};   // Taille suffisante pour "CH", "BH", etc.
    char *criteria = criteria_buf;
    char *N_especes = argv[1];
    const char *K_real = argv[2];
    char *percent = argv[3];

    int Strouve[treeAmount];

    /*for(int linej=0;linej<treeAmount;linej++){
        Strouve[linej]= 0;
    }*/

    copy(ctx.Strouve.begin(), ctx.Strouve.end(), Strouve);

    double *CHr, *Wr;
    CHr = new double [k_capacity+1];
    Wr = new double [k_capacity+1];

    //Ces variables sont utilisées pour déterminer combien de fois une boucle sera parcouru.
    double SSE=0,SSEref=0;

    int listr[k_capacity+1][nmax+1];

    for (int i=0; i<=k_capacity; i++){
        for (int j=0; j<=nmax; j++){
            listr[i][j] = 1;
        }
    }

    int *list,*no;        //list(nmax),no(nmax)
    list = new int [nmax+1];
    no = new int [nmax+1];

    for (int i=0; i<=nmax; i++){
        list[i] = 0;
        no[i] = 0;
    }

    int *howmany;        //howmany(k_capacity),
    howmany = new int [k_capacity+1];

    for (int i=0; i<=k_capacity; i++){
        howmany[i] = 0;
    }

//***********************  Read data file  **********************************

    setup_k_bounds(k_max, k_capacity, k_min, k1, k2, treeAmount, isBH, CHr, Wr);

    //--Read the data from files
    ReadData1(treeAmount,nmax,pmax);

    for(int i1=0; i1<treeAmount; i1++){
        for(int i2=0; i2<treeAmount; i2++){
            mat[i1][i2] = arrondir(mat[i1][i2],ROUNDING_PRECISION);
        }
    }

    double CH_new = MIN_CH_VALUE;
    double W_new = MAX_W_VALUE;

    int Sref [treeAmount];
    for (int j=0; j<treeAmount; j++){
        Sref[j]=0;
    }

    int nbInit =0;
    int nbFin =0;
    map <int, int> CH_conversion;
    map <int, int> W_conversion;

    //Est-ce que on as vraiment besoin que ça soit 3 variables séparées?
    int realk = 0;
    int CHk = 0, wk = 0;

    for (int i=0; i<=k_capacity; i++){
        howmany[i] = 0;
    }

    //iran is the number of random position
    for (int iran=1;iran<=nran; iran++) {
        CH_new = MIN_CH_VALUE;
        CHk = 0;
        wk = 0;
        realk = 0;

        if(iassign!=4){
            Assign(iran,treeAmount,k1,list,howmany,no,iassign,random_number);
        }
        // Big loop on number of groups, downwards from k1 to k2 (k1>=k2) - - - - - -

        //initialisation de Strouve de la liste realiser aleatoirement
        for (currentK=k1;currentK>=k2;currentK--){
            SSEref=INITIAL_SSE_REF;
            FO_new = MAX_FO_VALUE;
            W_new = MAX_FO_VALUE;

            //Je ne comprends pas pourquoi on fait une boucle ici, il n'est peut-être pas nécéssaire de faire ceci 100 fois.
            for (int nit=1;nit<=MAX_ITERATIONS;nit++){
                if(debug){
                    printf ("Iteration = %d",nit);
                    printf ("SSEref = %lf",SSEref);
                    for (int i=1;i<=treeAmount;i++){
                        printf ("%d",list[i]);
                    }
                    printf("\n");
                }

                // Compute distances to group centroids and assign objects to nearest one
                if(!isBH){
                    FO_new = FO_super_tree(treeAmount,k_capacity,mat,list,howmany,SSE,currentK);
                    CH_new = DistanceCH(treeAmount,k_capacity,mat,list,FO_new);
                    if(CH_new>CHr[currentK]){
                        CHr[currentK]=CH_new;
                        for (int i=1;i<=treeAmount;i++) {            //do 65 i=1,n
                            listr[currentK][i]=list[i];
                        }    //65    listr(currentK,i)=list(i)
                    }
                }else if(isBH){
                    FO_new = FO_W(treeAmount,k_capacity,mat,list,howmany,SSE,currentK);
                    W_new = DistanceW(treeAmount,k_capacity,list,FO_new);

                    if(W_new<Wr[currentK]){
                        Wr[currentK]=W_new;

                        for (int i=1;i<=treeAmount;i++){
                            listr[currentK][i]=list[i];
                        }
                    }
                }

                // Compute sum of squared error statistic (SSE) = within-group sum of squares

                if(fabs(SSEref-SSE)>(SSE/CONVERGENCE_THRESHOLD_DIVISOR)) {
                    SSEref=SSE;
                    if(nit==MAX_ITERATIONS) {
                        printf ("Convergence not reached in %d iterations.\n",MAX_ITERATIONS);
                    }
                }else{
                    nit = MAX_ITERATIONS + 1;   //Juste pour sortir d'un for loop
                }
            }

            // Compute the Calinski-Harabasz (1974) index 'CH' and

            // Concatenate the two closest groups before going to the next value of currentK
            int i1ref=1;        //i1ref=igr1
            int i2ref=currentK;        //i2ref=igr2

            //Group "i2ref" disappears
            for (int i=1;i<=treeAmount;i++){
                if(list[i]==i2ref){
                    list[i]=i1ref;
                }
                if(list[i]==i2ref) list[i]=list[i]-1;
            }        //70 continue

            howmany[i1ref]=howmany[i1ref]+howmany[i2ref];

            for (int k=(i2ref+1);k<=currentK;k++){
                howmany[k-1]=howmany[k];
            }

        }   //end for each k
        //--------------------------------------------------------------------
        // Affichage des organisation des groupes pour chaque nran (random start)
        //--------------------------------------------------------------------

        nbInit =0;
        nbFin =0;

        for(int i=0;i<tabIndices.size();i++){
            nbFin+=tabIndices.at(i);
            for (int j=nbInit; j<nbFin; j++){
                Sref[j]=i+1;
            }
            nbInit+=tabIndices.at(i);
        }

        if(!isBH){
            for (int k=k1;k>=k2;k--){
                if (CHr[k]>=CHr_max){
                    CHr_group=k;
                    CHr_max=CHr[k];

                    //Pour évider les clusters vides
                    realk=0;
                    CHk = 0;

                    for(int i=1; i<=k; i++){
                        for(int j=1; j<=treeAmount; j++){
                            if(listr[CHr_group][j]==i){
                                CHk++;
                                CH_conversion[i] = CHk;
                                realk++;
                                j = treeAmount + 1;
                            }
                        }
                    }

                    for(int iz=1; iz<=treeAmount; iz++){
                        listr[CHr_group][iz]=CH_conversion[listr[CHr_group][iz]];
                        Strouve[iz-1]=CH_conversion[listr[CHr_group][iz]];
                    }
                    CHr_group=realk;
                }
            }
        }else if(isBH){
            for (int k=k1;k>=k2;k--){
                if (Wr[k]<=W_min){
                    //pour connaitre le nombre de partition adéquate.
                    W_group=k;
                    W_min=Wr[k];

                    realk=0;
                    wk = 0;

                    for(int i=1; i<=k; i++){
                        for(int j=1; j<=treeAmount; j++){
                            if(listr[W_group][j]==i){
                                wk++;
                                W_conversion[i] = wk;
                                realk++;
                                j = treeAmount + 1;
                            }
                        }
                    }

                    for(int iz=1; iz<=treeAmount; iz++){
                        listr[W_group][iz]=W_conversion[listr[W_group][iz]];
                        Strouve[iz-1]=W_conversion[listr[W_group][iz]];
                    }
                    W_group=realk;
                }
            }
        }
    }  //fin random start

    // Print results
    if (!isBH){
        strcpy(criteria, "CH");
        conv2sameRef(Strouve,Sref,treeAmount);
        outStat(Strouve,Sref,criteria,treeAmount,N_especes,percent,K_real,CHr_group,CHr_max,monTableau);
    }else if(isBH){
        strcpy(criteria, "BH");
        conv2sameRef(Strouve,Sref,treeAmount);
        outStat(Strouve,Sref,criteria,treeAmount,N_especes,percent,K_real,W_group,W_min,monTableau);
    }

    // End timer
    const time_t tend2 = time(NULL);                // get the current calendar time

    // Compute execution time
    const double texec2 = difftime(tend2,tbegin2);    // tend-tbegin (result in second)
    fprintf (Output4,"%.3f;\n",texec2);

    // cleanup resources
    kmeans_cleanup(Output4, k_capacity, treeAmount, CHr, Wr,
        list, no, howmany);

    return 0;
}

int setup_k_bounds(int k_max, int k_capacity, int k_min, int &k1, int &k2, int treeAmount, bool isBH, double *&CHr, double *&Wr) {
    int warningCount = 0;
    int max_k1 = k_max;

    if (treeAmount<=k_min) max_k1=treeAmount-1;

    k1=max_k1;
    double facteur = 1.0;
    k2=k_min;
    if(!isBH){
        for (int i=0; i<=k_capacity; i++){
            CHr[i] = MIN_CH_VALUE;
        }
        if (k1<=2) {
            printf("*** Warning, not enough trees (k1:%d) k1 set to %d\n",k1,max_k1);
            k1=max_k1;
            warningCount++;
        }
        if (k_min<2){
            k2=2;
        }
    }else if(isBH){
        for (int i=0; i<=k_capacity; i++){
            Wr[i] = MAX_W_VALUE;
        }
        if (k_min<1){
            k2=1;
        }
        if (k1<=1) {
            printf("*** Warning, not enough trees (k1:%d) k1 set to %d\n",k1,max_k1);
            k1=max_k1;
            warningCount++;
        }
    }

    if (k1>k_capacity) {
        printf("*** Warning, limiting groups to %d \n",k_capacity);
        k1=max_k1-1;
        warningCount++;
    }

    return warningCount;
}

void kmeans_cleanup(FILE *Output4, int k_capacity, int treeAmount,
                    double *CHr, double *Wr,
                    int *list, int *no, int *howmany) {
    //Close output files
    if (Output4) fclose(Output4);

    delete [] CHr;
    delete [] Wr;
    delete [] list;
    delete [] no;
    delete [] howmany;
}

//      end
//************************End of Main

//******************************************************************************
//**********************************FUNCTIONS***********************************
//******************************************************************************

void ReadData1(int treeAmount,int nmax,int pmax){
    char *nameb;
    nameb = new char [MAX_FILENAME_LENGTH];

    if(treeAmount>nmax) {
        printf ("Too many objects. Use a sample of objects or recompile program to increase nmax.");                //     +'Too many objects. Use a sample of objects or recompile program.'
        exit(1);
    }

    if(treeAmount>pmax) {
        printf ("Too many variables. Use a sample of objects or recompile program to increase pmax.");                //     +'Too many objects. Use a sample of objects or recompile program.'
        exit(1);
    }

   //fclose(Input1);

    strcpy(nameb,"../output/stat.csv");
    if((Output4 = fopen(nameb,"a"))==NULL){
        printf("\n%s: result file open failed...",nameb);
        exit(1);
    }
    
    delete [] nameb;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

void Assign(int iran,int treeAmount,int k1,int* list,int* howmany,int* no,int iassign, int random_number){
    int ii=0, how=0, isum=0;
    char namea[MAX_PATH_LENGTH];

    if ((iassign==1) || (iassign==2)) {
        how=treeAmount/(k1*1.0);

        for (int k=1;k<=(k1-1);k++) {
            howmany[k]=how;
        }
        howmany[k1]=treeAmount-(k1-1)*how;

        for (int k=1;k<=k1;k++) {
            for (int currentK=1;currentK<=howmany[k];currentK++) {
               ii++;
               list[ii]=k;
            }
        }

        if(iassign!=1) {
            // Assign objects at random to the groups
            if(iran==1) {
                for (int i=1;i<=(random_number+100);i++) {
                    rand()/(1.0*(rand() % RAND_MAX_VALUE));
                }
            }
            Permute(treeAmount,list);
        }
    } else if (iassign==3) {
        // Read file of group assignments.
        // First line: how many objects in each group?
        // Then, read members of each group on one line (list of object numbers).

        printf ("Name of file of group assignments?");        //60 write(*,*) 'Name of file of group assignments?'
        scanf ("%s",namea);        //read(*,*) namea

        FILE *Input3;
        if ((Input3 = fopen(namea,"r"))==0) {
            printf("\n %s :Open Failed....",namea); exit(1);
        }
        printf ("File of group assignments: %s\n",namea);

        for (int k=1;k<=k1;k++) {
            fscanf(Input3,"%d",&howmany[k]);
        }

        isum=0;
        for (int k=1;k<=k1;k++) {
            isum=isum+howmany[k];
        }

        if(isum!=treeAmount) {
            printf("Objects assigned to groups do not sum to n.");
            exit(1);
        }

        for (int i=1;i<=treeAmount;i++) {
            list[i]=-1;
        }

        for (int k=1;k<=k1;k++){
            for (int i=1;i<=howmany[k];i++) {
                fscanf(Input3, "%d", &no[i]);
            }
            for (int i=1;i<=howmany[k];i++) {
                list[no[i]]=k;
            }
        }

        for (int i=1;i<=treeAmount;i++) {
            if(list[i]==-1) {
                printf("Overlapping assignments to groups.");
                exit(1);
            }
        }
        fclose(Input3);
    } else {
        printf("Wrong perameter <iassign> in function <Assign>.");
        exit(1);
    }
    return;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

//***********Permute

// This subroutine permutes the first 'n' values of vector 'iordre' at random
// in an equiprobable way. This property has been checked through intensive
// simulations.

void Permute(int n,int *iordre){
    // On parcourt le tableau de la dernière position vers la deuxième.
    // À chaque étape, un élément est échangé avec un élément choisi aléatoirement parmi les positions restantes.

    for (int m = n; m >= 2; --m) {
        // On genere un indice aleatoire j compris entre 1 et m inclus
        int j = 1 + (rand() % m);

        // On echange les éléments situés aux positions m et j
        int tmp = iordre[m];
        iordre[m] = iordre[j];
        iordre[j] = tmp;
    }
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

//compute rand index
double f_RI(int Strouve[],int Sref[],int treeAmount){
    double comb = 1.0;

    for (int i=treeAmount; i>=(treeAmount-2+1); i--) {
        comb*=i;
    }

    comb/=2.0;

    double a=0.0;
    double b=0.0;

    for (int i=0; i<treeAmount-1; i++){
        for (int j=i+1; j<treeAmount; j++){
            if(Sref[i]!=Sref[j]){
                if(Strouve[i]!=Strouve[j]){
                    b++;
                }
            }else{
                if(Strouve[i]==Strouve[j]){
                    a++;
                }
            }
        }
    }
    double RI = (a+b)/(comb*1.0);

    return RI;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

//compute adjusted rand index
double f_ARI(int Strouve[],int Sref[],const char *K_real,int group,int treeAmount){
    int kReal = atoi(K_real);
    if(kReal <= 0){
        // If K_real is not a valid positive integer (e.g. "?"),
        // fall back to using the detected 'group' value to avoid
        // creating zero-sized VLA and prevent out-of-bounds access.
        kReal = group;
    }
    int tabCongruence [kReal+1][group+1];
    int sumLigne [kReal+1];
    int sumColonne[group+1];

    //initialisation du tableau sumLigne à 0
    for(int i=0;i<=kReal;i++){
        sumLigne[i]=0;
    }

    //initialisation du tableau sumColonne à 0
    for(int i=0;i<=group;i++){
        sumColonne[i]=0;
    }

    //initialisation du tableau des congruence à 0
    for(int i=0;i<=kReal;i++){
        for(int j=0;j<=group;j++){
            tabCongruence[i][j]=0;
        }
    }


    for(int i=0;i<treeAmount;i++){
        tabCongruence[Sref[i]][Strouve[i]]++;
    }

    for(int i=1;i<=kReal;i++){
        for(int j=1;j<=group;j++){
            sumLigne[i]+=tabCongruence[i][j];
            sumColonne[j]+=tabCongruence[i][j];
        }
    }

    double a=0.0;
    double b=0.0;
    double c=0.0;

    double comb = 1.0;

    for (int i=treeAmount; i>=(treeAmount-2+1); i--) {
        comb*=i;
    }

    comb/=2.0;

    for (int i=0; i<treeAmount-1; i++){
        for (int j=i+1; j<treeAmount; j++){
            if(Sref[i]!=Sref[j]){
                if(Strouve[i]==Strouve[j]){
                    c++;
                }
            }else{
                if(Strouve[i]==Strouve[j]){
                    a++;
                }else{
                    b++;
                }
            }
        }
    }

    double ARI = 0.0;

    if(a*2.0==((b+a)+(c+a))){
        ARI = 1.0;
    }else{
        ARI = a - ((b+a)*(c+a))/(comb*1.0);
        ARI=ARI/((((b+a)+(c+a))/2.0)-(((b+a)*(c+a))/(comb*1.0)));
    }

    return ARI;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

// Modification Centroids: this whole subroutine

//stat output
void outStat(int Strouve[],int Sref[],char *criteria,int treeAmount,char *N_especes,char *percent,const char *K_real,int group,double score, vector <string> monTableau){
    //Compute Rand index between Strouve and Sref
    double RI = f_RI(Strouve,Sref,treeAmount);

    //Compute Rand index Adjusted between Strouve and Sref
    double ARI = f_ARI(Strouve,Sref,K_real,group,treeAmount);

    fprintf (Output4,"%s;",criteria);
    fprintf (Output4,"%i;",treeAmount);
    fprintf (Output4,"%s;",N_especes);
    fprintf (Output4,"%s;",percent);
    fprintf (Output4,"%s;",K_real);
    fprintf (Output4,"%d;",group);
    int diff = atoi(K_real)-group;
    diff = fabs(diff);
    const int const_K_real = atoi(K_real);
    const int const_group = group;

    double max_k = max(const_K_real,const_group);
    double diff_norm = (diff*1.0)/(max_k*1.0);
    fprintf (Output4,"%i;",diff);
    fprintf (Output4,"%.3f;",diff_norm);
    if(atoi(K_real) == group){
        fprintf (Output4,"%i;",1);
    }else{
        fprintf (Output4,"%i;",0);
    }
    fprintf (Output4,"%.3f;",RI);
    fprintf (Output4,"%.3f;",ARI);
    fprintf (Output4,"%.3f;",score);

    fprintf (Output4,"part(");
    for (int p=1; p<=treeAmount; p++){
        if(p==treeAmount){
            fprintf (Output4,"%i%s",Strouve[p-1]," ");
        }else{
            fprintf (Output4,"%i%s",Strouve[p-1]," <> ");
        }

    }
    fprintf (Output4,");");
    
    //output.txt file
    //composition of each cluster and each element
    ofstream myfile;
    myfile.open ("../output/output.txt");
    myfile << criteria;
    myfile << "\n";
    for(int n_cl=1; n_cl<=group; n_cl++){
        myfile << "Cluster # ";
        myfile << n_cl;
        myfile << "\n";
        myfile << "Cluster content:";
        myfile << "\n";
        for (int p=1; p<=treeAmount; p++){
            if(Strouve[p-1]==n_cl){
                myfile << "\tTree #\t: T";
                myfile << p;
                myfile << "\t";
                myfile << monTableau[p-1];
                myfile << "\n";
            }
        }
        myfile<<"\n";
    }
    myfile<<"\n";
    myfile.close();
    
    std::cout<<"1) stat.csv - for clustering statistics;"<< std::endl;
    std::cout<<"2) output.txt - for cluster content."<< std::endl;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double FO_super_tree(int treeAmount, int k_capacity, double** mat, int* list, int* howmany, double &SSE, int currentK){
    // clusterK_same[k] stocke la somme des distances RF internes (ou vers un représentant)
    // utilisée pour calculer la contribution du cluster k à la fonction objectif.
    double *clusterK_same = new double[k_capacity + 1];

    // nk_CH[k] stocke le nombre d'arbres actuellement assignés au cluster k.
    int *nk_CH = new int[k_capacity + 1];

    SSE = 0.0;

    for (int k = 1; k <= k_capacity; ++k) {
        nk_CH[k] = 0;
        clusterK_same[k] = 0.0;
        howmany[k] = 0;
    }

    //boucle sur les arbres
    // Compte le nombre d'arbres dans chaque cluster à partir des affectations list[1..n].
    for (int i = 1; i <= treeAmount; ++i) {
        int g = list[i];
        if (g >= 1 && g <= k_capacity) {
            nk_CH[g]++;
        }
    }

    // Met à jour howmany[k] avec les effectifs des clusters (1..currentK).
    for (int k = 1; k <= currentK; ++k) {
        howmany[k] = nk_CH[k];
    }

    // ------------------------------------------------------------
    // Calcul de clusterK_same[k]
    // Si withConsensus == true : on utilise un "médoïde" :
    //   - pour chaque cluster k, on choisit l'arbre du cluster qui minimise
    //     la somme des distances RF vers tous les autres arbres du cluster.
    //   - clusterK_same[k] devient cette somme minimale.
    //
    // Sinon : on calcule la somme des distances RF sur toutes les paires
    // d'arbres du cluster (intra-cluster).
    // ------------------------------------------------------------
    if (withConsensus) {
        for (int k = 1; k <= currentK; ++k) {
            // Si le cluster contient 0 ou 1 arbre, sa contribution est nulle.
            if (nk_CH[k] <= 1) {
                clusterK_same[k] = 0.0;
                continue;
            } else {
            double bestSum = 1e100;

                // On teste chaque arbre i du cluster comme représentant candidat.
                for (int i = 1; i <= treeAmount; ++i) {
                    if (list[i] == k) {
                        double sumDist = 0.0;

                        // Somme des distances RF entre i et tous les autres arbres j du même cluster.
                        for (int j = 1; j <= treeAmount; ++j) {
                            if ((i != j)&&(list[j] == k)) {
                                sumDist += mat[i - 1][j - 1];
                            }
                        }

                        // On conserve le candidat qui minimise la somme des distances.
                        if (sumDist < bestSum) {
                            bestSum = sumDist;
                        }
                    }
                }
                clusterK_same[k] = bestSum;
            }
        }
    } else {
        //boucle sur les arbres
        // Somme des distances RF sur toutes les paires (i,j) dans le même cluster.
        for (int i = 1; i < treeAmount; ++i) {
            int cluster_i = list[i];
            for (int j = i + 1; j <= treeAmount; ++j) {
                if (list[j] == cluster_i) {
                    clusterK_same[cluster_i] += mat[i - 1][j - 1];
                }
            }
        }
    }

    // ------------------------------------------------------------
    // Calcul de la fonction objectif FO_old :
    // Pour chaque cluster k, on ajoute une contribution normalisée par le nombre d'arbres.
    // ------------------------------------------------------------
    double FO_old = 0.0;
    for (int k = 1; k <= currentK; ++k) {
        if (nk_CH[k] > 1) {
            FO_old += (clusterK_same[k] / (1.0 * nk_CH[k]));
        }
    }

    // Par défaut, on retourne la valeur actuelle de la fonction objectif.
    // (Ton code original tente ensuite d'améliorer l'affectation en déplaçant des arbres.)
    double Dref = FO_old;

    delete [] clusterK_same;
    delete [] nk_CH;

    return Dref;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double DistanceCH(int treeAmount,int k_capacity,double** mat,int* list,double FO_new){
    double SSB = 0.0;
    double SSW = 0.0;
    double dist_all = 0.0;
    double RF;
    double distance_total = 0.0;
    int *nk_CH = new int [k_capacity+1];
    int k_cluster = 0;

    for(int k=1;k<=k_capacity; k++){
        nk_CH[k]=0;
    }

    // On parcourt tous les arbres (1..n) pour compter combien appartiennent à chaque cluster.
    for (int i = 1; i <= treeAmount; ++i) {
        int g = list[i];            // g = numéro de cluster attribué à l'objet i
        if (g >= 1 && g <= k_capacity) {  // on vérifie que g est dans les bornes du tableau nk_W
            nk_CH[g]++;             // on incrémente le compteur du cluster g
        }
    }

    for(int k=1;k<=k_capacity; k++){
        if(nk_CH[k]!=0){
            k_cluster++;
        }
    }

    //compute dist_all
    for (int i=1;i<treeAmount;i++){
        for (int j=i+1;j<=treeAmount;j++){
            RF = mat[i-1][j-1];
            dist_all += RF;
        }
    }

    dist_all = dist_all/(1.0*treeAmount);

    //compute SSW
    SSW = FO_new;

    //compute SSB
    SSB = dist_all - SSW;

    if((fabs(SSW)>0.000001) && (k_cluster>1)){
        distance_total=(SSB/(1.0*SSW))*((1.0*treeAmount-k_cluster)/((1.0*k_cluster)-1.0));
    }
    if(fabs(SSW)<=0.000001  && (k_cluster>1)){
        distance_total=10000000.0*SSB*((treeAmount-k_cluster)/((1.0*k_cluster)-1.0));
    }

    delete [] nk_CH;

    return distance_total;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double FO_W(int treeAmount,int k_capacity,double** mat,int* list,int* howmany,double &SSE,int currentK){
    double *clusterK_same = new double [k_capacity+1];
    int *nk_W = new int [k_capacity+1];
    int cluster_k = 0;
    double RF = 0.0;
    double Dref = 0;       //Real*8 Dref,D1,SSE,weight(pmax)
    int    kref = 0;        //Integer list(nmax),howmany(k_capacity),kref
    // Compute squared distances to group centroids. Assign objects to nearest one
    SSE=0;

    int k_source = 0;
    int nb_cluster_dest = 0;
    int nb_cluster_source = 0;
    double FO_old = 0.0;
    double FO_new = 0.0;
    double tmp_calc_dest = 0.0;
    double tmp_calc_source = 0.0;

    for(int k=1;k<=k_capacity; k++){
        nk_W[k]=0;
        clusterK_same[k]=0.0;
    }

    // On parcourt tous les arbres (1..n) pour compter combien appartiennent à chaque cluster.
    // Comptage sécurisé : évite nk_W[list[i]] hors bornes si list[i] est invalide.
    for (int i = 1; i <= treeAmount; ++i) {
        int g = list[i];            // g = numéro de cluster attribué à l'objet i
        if (g >= 1 && g <= k_capacity) {  // on vérifie que g est dans les bornes du tableau nk_W
            nk_W[g]++;              // on incrémente le compteur du cluster g
        }
    }

    //calcul distance
    //compute for each cluster initially, SSW value (intra groupe distance)
    //compute SSW
    for (int i=1;i<treeAmount;i++){
        cluster_k=list[i];
        for (int j=i+1;j<=treeAmount;j++){
            if (list[j]==cluster_k){
                RF = mat[i-1][j-1];
                clusterK_same[cluster_k]+=RF;
            }
        }
    }

    for (int k=1;k<=currentK;k++){
        if(nk_W[k]>1){
            FO_old += ((2.0*clusterK_same[k])/(1.0*nk_W[k]*(nk_W[k]-1)));
        }
    }

    if(currentK==1){
        Dref=FO_old;
    }else{
        //boucle sur les arbres
        for (int i=1;i<=treeAmount; i++){
            if(nk_W[list[i]]>1){
                for (int k=1;k<=currentK;k++){
                    //Calcul de la distance RF de chaque point i
                    // et assignation du point i au bon cluster
                    // Compute a RF distance to the centroid k

                    //test si le point i n'appartenait pas initiallement à k
                    //k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
                    if((list[i]!=k) && (nk_W[list[i]]>1)){
                        //Pour le cluster source
                        k_source = list[i];
                        nb_cluster_source = nk_W[k_source];

                        tmp_calc_source = clusterK_same[k_source];
                        nb_cluster_source -=1;
                        if (nb_cluster_source==1){
                            tmp_calc_source = 0.0;
                        }else{
                            for(int j=1;j<=treeAmount; j++){
                                if(list[j]==k_source && i!=j){
                                    tmp_calc_source -= mat[i-1][j-1];
                                }
                            }
                        }

                        nb_cluster_dest = nk_W[k];

                        if((nk_W[k]>1) && (nk_W[list[i]]>1)){
                            FO_new = FO_old - (2.0*clusterK_same[k]/(1.0*nk_W[k]*(nk_W[k]-1)));
                            FO_new = FO_new - (2.0*clusterK_same[list[i]]/(1.0*nk_W[list[i]]*(nk_W[list[i]]-1)));
                        }
                        else if(nk_W[list[i]]>1){
                            FO_new = FO_old - (2.0*clusterK_same[list[i]]/(1.0*nk_W[list[i]]*(nk_W[list[i]]-1)));
                        }
                        else if(nk_W[k]>1){
                            FO_new = FO_old - (2.0*clusterK_same[k]/(1.0*nk_W[k]*(nk_W[k]-1)));
                        }
                        else{
                            FO_new = FO_old;
                        }
                        tmp_calc_dest = clusterK_same[k];
                        for(int j=1;j<=treeAmount; j++){
                            if(list[j]==k){
                                tmp_calc_dest += mat[i-1][j-1];
                            }
                        }
                        nb_cluster_dest +=1;
                        if(nb_cluster_dest>=1){
                            FO_new = FO_new + (2.0*tmp_calc_dest/(1.0*nb_cluster_dest*(nb_cluster_dest-1)));
                        }

                        //Pour le cluster source
                        if(nb_cluster_source>=1){
                            FO_new = FO_new + (2.0*tmp_calc_source/(1.0*nb_cluster_source*(nb_cluster_source-1)));
                        }

                        if(FO_new<FO_old){
                            Dref=FO_new;
                            kref=k;

                            //mise à jour de nk_W[]
                            nk_W[k] = nb_cluster_dest;
                            nk_W[list[i]] = nb_cluster_source;

                            //mise à jour de la distance intra-groupe des deux clusters modifiés
                            clusterK_same[list[i]] = tmp_calc_source;
                            clusterK_same[k] = tmp_calc_dest;

                            howmany[k] = nb_cluster_dest;
                            howmany[list[i]] = nb_cluster_source;

                            //calcul SSE
                            SSE=SSE+Dref;

                            //mise à jour de la fonction objective FO_old
                            FO_old = FO_new;

                            //déplacement de l'arbre
                            //mise à jour la liste de distribution des elements
                            list[i] = k;
                        }else{
                            Dref=FO_old;
                        }
                    }
                }
            }
        }
    }

    delete [] clusterK_same;
    delete [] nk_W;

    return Dref;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double DistanceW(int treeAmount, int k_capacity, int* list, double FO_new){
    double distance_total = 100000000.0;
    double *clusterK_same = new double [k_capacity+1];
    int *nk_W = new int [k_capacity+1];
    int k_cluster = 0;

    for(int k=1;k<=k_capacity; k++){
        nk_W[k]=0;
        clusterK_same[k]=0.0;
    }

    // Comptage sécurisé : évite nk_W[list[i]] hors bornes si list[i] est invalide.
    for (int i = 1; i <= treeAmount; ++i) {
        int g = list[i];            // g = numéro de cluster attribué à l'objet i
        if (g >= 1 && g <= k_capacity) {  // on vérifie que g est dans les bornes du tableau nk_W
            nk_W[g]++;              // on incrémente le compteur du cluster g
        }
    }

    for(int k=1;k<=k_capacity; k++){
        if(nk_W[k]!=0){
            k_cluster++;
        }
    }

    if(k_cluster!=treeAmount){
        // distance_total=(FO_new/(1.0*(n-k_cluster)));
        if(k_cluster!=0) {
            distance_total=(FO_new/k_cluster);
        }
    }

    delete [] clusterK_same;
    delete [] nk_W;

    return distance_total;

}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double arrondir(double num,int digits){
    return floor(num*pow(10,digits)+0.5)/(1.0*pow(10,digits));
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

//Convert the partition found by the same number cluster that the partition ref
void conv2sameRef(int *Strouve,int *Sref, int treeAmount){
    int k = 0;

    std::cout<<"Number of trees in the input file: "<<treeAmount<< std::endl;
    std::cout<<"Partition found: "<< std::endl;
    std::cout<<Strouve[0];
    for(int i=1; i<treeAmount; i++){
        std::cout<<" <> "<<Strouve[i];
    }
    std::cout<<std::endl;

    //To know the number of cluster
    for(int i=0; i<treeAmount; i++){
        if(Strouve[i]>k){
            k=Strouve[i];
        }
    }

    std::cout << "Number of clusters (K) found: " << k << std::endl;

    // On utilise k+1 cases car les clusters sont numérotés de 1 à k (l'indice 0 n'est pas utilisé).
    int *nk_trouve = new int[k + 1];
    int *nk_ref    = new int[k + 1];

    // On met à 0 le nombre d'arbres dans chaque cluster (1..k).
    for (int c = 0; c <= k; c++) {
        nk_trouve[c] = 0;
        nk_ref[c] = 0;
    }

    // Pour chaque arbre (0..treeAmount-1), on incrémente le compteur du cluster auquel il appartient.
    for (int i = 0; i < treeAmount; i++) {
        nk_trouve[Strouve[i]]++;
        nk_ref[Sref[i]]++;
    }

    // Ici, nk_trouve[c] = nombre d'arbres assignés au cluster c dans la partition trouvée.
    // nk_ref[c]    = nombre d'arbres assignés au cluster c dans la partition de référence.

    delete [] nk_trouve;
    delete [] nk_ref;
}