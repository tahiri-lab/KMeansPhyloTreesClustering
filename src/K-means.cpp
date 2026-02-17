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
const double INITIAL_MIN_W = 1e9;               // Initial minimum W value for optimization
const double INITIAL_MAX_CH = -1e10;            // Initial maximum CH value
const int MAX_ITERATIONS = 100;                 // Maximum K-means iterations per cluster size
const int RAND_MAX_VALUE = 32767;               // Maximum random number from rand()
const int MAX_FILENAME_LENGTH = 300;            // Maximum length for filename strings
const int MAX_PATH_LENGTH = 255;                // Maximum length for file paths
const int DISTANCE_ARRAY_SIZE = 4;              // Size of Robinson-Foulds distance array
const int ROUNDING_PRECISION = 3;               // Decimal places for rounding distance values
const int CONVERGENCE_THRESHOLD_DIVISOR = 1000; // Divisor for convergence threshold check
const double MIN_DISTANCE = 1000000.0;          // Minimum distance value for clustering
// =============================================================================================================

FILE *Output4;


//  Parameters of the program:
//  n    = number of observations (ex. number of trees)
//  nmax = maximum number of observations (default: 10000)
//  p    = number of variables     (ex. number of variables describing each trees)
//  pmax = maximum number of variables  (default: 10000)
//  k    = number of groups (centroids)
//  kmax = maximum number of groups
//  kk   = ???
//  niter = maximum iteration for convergeance of centroid (fixed=100)
//  Parameter (nmax=100000,pmax=10,kmax=10)
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
//   ishort[p] = vector containing the position(p) of valide variables (containing non zero value for weight -- see ReadData).
//   nobest
//   var       = double [kmax+1]
//   vect      = double [pmax+1]
//   mean[p]      = vector of mean value not weigthed for the variable [p] double [pmax+1]
//   Dvec      = double [kmax+1]
//   CHr       = vector[k] of the best CH this k
//   BHr       = vector[k] of the best BH this k
//   Silr      = vector[k] of the best Silr this k
//   LogSSr    = vector[k] of the best LogSS this k
//   Wr        = vector[k] of the best W this k
//   SSEr      = double [kmax+1]  sum of squared error statistic (SSE) = within-group sum of squares
//   listr     = int[kmax+1][nmax+1];
//   howmanyr  = int[kmax+1][nmax+1];
// sx = new double*[kmax+1];
//  sx2 = new double*[kmax+1];

//   no = new int [nmax+1];
//    iordre = new int [nmax+1];
//   nobest = new int [kmax+1];
//   nnitr = new int [kmax+1];
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

int main_kmeans(char **argv,vector <string> monTableau, double ** mat, double ** n_identique,double ** Ww, vector<int> tabIndices, int intParam, int *n_leaves,int k_min, int k_max, int alpha)
{
    //*****************Define variables******************************************//
    // Variables
    map<int,string>  mapIndicesTreesFinal;
    vector <string> indicesTrees;
    time_t tbegin2,tend2;
        double texec2=0.;

    double W = 0.0;
    double CH = MIN_CH_VALUE;
    bool use_weight=false;

    double CHr_max=INITIAL_MAX_CH;
    int CHr_group=0;

    double W_min=MAX_W_VALUE;
    double W_max=MIN_CH_VALUE;
    int W_group=0;
    double FO_new = MAX_FO_VALUE;

        // Start timer
        tbegin2=time(NULL);                // get the current calendar time

    int N = int (monTableau.size()); //quantity of initial tree
    int i=0, j=0;        //Counters
    int n=N, p=N;        //,pmax,kmax; //      Integer p,pmax,kmax
    int ntran=0, iseed=0, niter=0, kk=0, nit=0;        //added declarations for variables
    int nnit=0, k=0, i1ref=0, i2ref=0;        //added declarations for variables
    int idebug=0 ; // 0, no debug, 1 debug
    int k1=0, k2=0;  //added declarations for variables
    int hard_max_k=0; //--Setting the max k1

        int random_number=100; //--Fixed random number
    int iassign=2;  // 1 equal, 2 random
        int iran=100;   //--Number of random position
        int nran=100;  //--Number of Random start VM

    int nmax=N;    //--Maximum number of object -Parameter (nmax=10000,pmax=250,kmax=100)
    int pmax=N;      //--Maximum data point (variable))
     int kmax=N;      // Maximum number of groups

    char *criteria = argv[0];
    char *N_especes = argv[1];
    const char *K_real = argv[2];
    char *percent = argv[3];

    int Strouve[N];

    for(int linej=0;linej<N;linej++){
        Strouve[linej]= 0;
    }

    double    **sx,**sx2,**xbar,**var;    //sx(kmax,pmax),sx2(kmax,pmax),xbar(kmax,pmax),var(kmax,pmax)
    double **tree_cluster_leaves = new double *[n];
    for(int i=0;i<n;i++){
        tree_cluster_leaves[i]=new double [DISTANCE_ARRAY_SIZE];
    }

    sx = new double*[kmax+1];
    sx2 = new double*[kmax+1];
    xbar = new double*[kmax+1];
    var = new double*[kmax+1];

    for (i=0;i<=kmax;i++){
        sx[i] = new double[pmax+1];
        sx2[i] = new double[pmax+1];
        xbar[i] = new double[pmax+1];
        var[i] = new double[pmax+1];
    }

    //variables centroids (trees and indices)
    vector <string> centroid_k;
    vector <int> centroid_k_pos;
    string centroid_C_min = "";

    double *distances_RF_norm = new double[DISTANCE_ARRAY_SIZE];

    for(int linej=0;linej<DISTANCE_ARRAY_SIZE;linej++){
        distances_RF_norm[linej]= 0.0;
    }

    for (i=0; i<=kmax; i++){
        for (j=0; j<=pmax; j++){
            sx[i][j] = 0.0;
            sx2[i][j] = 0.0;
            xbar[i][j] = 0.0;
            var[i][j] = 0.0;
        }
    }


    // double *Dvec,*CHr, *BHr,*SSEr, *Silr, *LogSSr, *Wr, *diff_W, *V_W, *Wr_ln, *Gapr;
    double *Dvec,*SSEr,*diff_W, *V_W, *Wr_ln, *CHr, *Wr;
    Dvec = new double [kmax+1];
    SSEr = new double [kmax+1];

    CHr = new double [kmax+1];
    Wr = new double [kmax+1];
    Wr_ln = new double [kmax+1];
    diff_W = new double [kmax+1];
    V_W = new double [kmax+1];


    for (i=0; i<=kmax; i++){
        Dvec[i] = 0.0;
        SSEr[i] = 0.0;
    }

    double *vect,*mean,*weight;        //vect(pmax),mean(pmax),weight(pmax),
    vect = new double [pmax+1];
    mean = new double [pmax+1];
    weight = new double [pmax+1];
    for (i=0; i<=pmax; i++){
        vect[i] = 0.0;
        mean[i] = 0.0;
        weight[i] = 0.0;
    }

    double D1=0,Dref=0,SSE=0,SSEref=0,SST=0;

    int **listr;                    //listr(kmax,nmax),
    listr = new int*[kmax+1];
    for (i=0;i<=kmax;i++){
        listr[i] = new int [nmax+1];
    }

    for (i=0; i<=kmax; i++){
        for (j=0; j<=nmax; j++){
            listr[i][j] = 1;
        }
    }

    int    **howmanyr;        //howmanyr(kmax,kmax)
    howmanyr = new int*[kmax+1];
    for (i=0;i<=kmax;i++){
        howmanyr[i] = new int [kmax+1];
    }

    for (i=0; i<=kmax; i++){
        for (j=0; j<=kmax; j++){
            howmanyr[i][j] = 0;
        }
    }


    int *list,*no, *iordre;        //list(nmax),no(nmax), iordre(nmax)
    list = new int [nmax+1];
    no = new int [nmax+1];
    iordre = new int [nmax+1];

    for (i=0; i<=nmax; i++){
        list[i] = 0;
        no[i] = 0;
        iordre[i] = 0;
    }

    int *howmany,*nobest/*, *nobestSilhouette, *nobestLogSS ,*nobestCH, *nobestBH, *nobestW */, *nnitr;        //howmany(kmax),,nobest(kmax), nnitr(kmax);
    howmany = new int [kmax+1];
    nobest = new int [kmax+1];
    nnitr = new int [kmax+1];

    for (i=0; i<=kmax; i++){
        howmany[i] = 0;
        nobest[i] = 0;
        nnitr[i] = 0;
    }

    int *ishort;            //ishort(pmax);
    ishort = new int [pmax+1];

    for (i=0; i<=pmax; i++){
        ishort[i] = 0;
    }


//  Modification Centroids: add ",nameb" to next line
    char *nameb;
    nameb = new char [MAX_FILENAME_LENGTH];


//***********************  Read data file  **********************************

       int max_k1 = k_max;

    if (N<=k_min) max_k1=N-1;

    k1=max_k1;
    double facteur = 1.0;
    k2=k_min;
    if(intParam==1){
        for (i=0; i<=kmax; i++){
            CHr[i] = MIN_CH_VALUE;
        }

        if (k1<=2) {
            printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
            k1=max_k1;
        }
        if (k_min<2){
            k2=2;
        }
    }else if(intParam==2){
        for (i=0; i<=kmax; i++){
            Wr[i] = MAX_W_VALUE;
            Wr_ln[i] = MIN_CH_VALUE;
            diff_W[i] = 0.0;
            V_W[i] = 0.0;
        }
        if (k_min<1){
            k2=1;
        }
        if (k1<=1) {
            printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
            k1=max_k1;
        }
    }

    if (k1>kmax) {
      printf("*** Warning, limiting groups to %d \n",kmax);
      k1=max_k1-1;
    }

    if (hard_max_k!=0) {
      k1=max_k1;
    }

     //--Read the data from files
    ReadData1(n,nmax,p,pmax,mat,ishort,weight,mean,ntran,nameb,N);    //Call ReadData11(n,nmax,p,pmax,mat,coord,ishort,w,mean,ntran,namea)

    CompSST(n,nmax,p,pmax,mat,weight,ishort,SST);

    for(int i1=0; i1<n; i1++){
        for(int i2=0; i2<n; i2++){
            mat[i1][i2] = arrondir(mat[i1][i2],ROUNDING_PRECISION);
        }
    }

    // Compute vector 'mean' of overall means
    for (j=1;j<=p;j++){
        mean[j]=0;
    }

    for (i=1;i<=n;i++){
        for (j=1;j<=p;j++){
            mean[j]=mean[j]+mat[i-1][ishort[j]-1];
        }
    }

    for (j=1;j<=p;j++){
        mean[j]=mean[j]/(n*1.0);//18 mean(j)=mean(j)/dfloat(n)
    }

    iseed=0;
    double CH_new = MIN_CH_VALUE;
    double W_new = MAX_W_VALUE;

    //initialize Sref
    int Sref [N];
    for (int j=0; j<N; j++){
        Sref[j]=0;
    }

    int number_cluster = 0;

    int nbInit =0;
    int nbFin =0;
    map <int, int> CH_conversion;
    map <int, int> W_conversion;

    int realk = 0;
    int unique = 0;
    int CHk = 0;
    int wk = 0;
    for (i=0; i<=kmax; i++){
        howmany[i] = 0;
        nobest[i] = 0;
        nnitr[i] = 0;
    }

    int *nk = new int [kmax+1];

    //initialisation des variables
    for(int k=1;k<=kmax; k++){
        nk[k]=0;
    }

    for (iran=1;iran<=nran; iran++) {
        CH_new = MIN_CH_VALUE;
        CHk = 0;
        use_weight = false;
        wk = 0;
        realk = 0;
        unique = 0;
        number_cluster = 0;


        if(iassign!=4){
            Assign(iran,n,nmax,k1,kmax,list,howmany,no,idebug,iassign,iseed, random_number);
        }
        // Big loop on number of groups, downwards from k1 to k2 (k1>=k2) - - - - - -
        niter=MAX_ITERATIONS; //changed VM

        //initialisation de Strouve de la liste realiser aleatoirement
        for (kk=k1;kk>=k2;kk--){
            SSEref=INITIAL_SSE_REF;
            W = MAX_W_VALUE;
            CH = MIN_CH_VALUE;
            FO_new = MAX_FO_VALUE;
            W_new = MAX_FO_VALUE;

            for (nit=1;nit<=niter;nit++){
                if(idebug==1){
                    printf ("Iteration = %d",nit);
                    printf ("SSEref = %lf",SSEref);
                    for (i=1;i<=n;i++){
                        printf ("%d",list[i]);
                    }
                }
                nnit=nit;

                // Compute distances to group centroids and assign objects to nearest one
                if(intParam==1){
                    FO_new = FO_super_tree(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau, alpha);
                }else if(intParam==2){
                    FO_new = FO_W(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau);
                }

                number_cluster = 0;

                if(intParam==1){
                    CH_new = DistanceCH(n,kmax,mat,list,Ww,FO_new,facteur);
                    if(CH_new>CHr[kk]){
                        SSEr[kk]=SSE;        //SSEr(kk)=SSE
                        nobest[kk]=iran;    //nobest(kk)=iran

                        nnitr[kk]=nnit;    //nnitr(kk)=nnit
                        CH=CH_new;
                        CHr[kk]=CH;  //CHr(kk)=CH
                        for (int i=1;i<=n;i++)            //do 65 i=1,n
                        {
                            listr[kk][i]=list[i];
                        }    //65    listr(kk,i)=list(i)

                        for (int i=1;i<=kk;i++)                //do 67 i=1,kk
                          {howmanyr[kk][i]=howmany[i];}    //67    howmanyr(kk,i)=howmany(i)
                    }
                }else if(intParam==2){
                    W_new = DistanceW(n,kmax,mat,list,Ww,FO_new,facteur);

                    if(W_new<Wr[kk]){
                        SSEr[kk]=SSE;        //SSEr(kk)=SSE
                        nobest[kk]=iran;    //nobest(kk)=iran
                        nnitr[kk]=nnit;    //nnitr(kk)=nnit

                        W=W_new;
                        Wr[kk]=W;

                        for (int i=1;i<=n;i++){
                            listr[kk][i]=list[i];
                        }

                        for (int i=1;i<=kk;i++){
                            howmanyr[kk][i]=howmany[i];
                        }
                    }
                }

                  // Compute sum of squared error statistic (SSE) = within-group sum of squares

                if(fabs(SSEref-SSE)>(SSE/CONVERGENCE_THRESHOLD_DIVISOR))            //if(dabs(SSEref-SSE).gt.SSE/1000.0) then
                {
                    SSEref=SSE;//SSEref=SSE
                }else{
                    goto m60;
                }

            }

            // Compute the Calinski-Harabasz (1974) index 'CH' and
            /* printf ("Convergence not reached in %d iterations.",niter);// write(*,*) 'Convergence not reached in ',niter,' iterations.' */
m60:
            // Concatenate the two closest groups before going to the next value of kk
            Dref=MIN_DISTANCE;    //Dref=1000000.0
            D1=0.0;        //D1=0.0
            i1ref=1;        //i1ref=igr1
            i2ref=kk;        //i2ref=igr2

            //Group "i2ref" disappears
            for (i=1;i<=n;i++){
                if(list[i]==i2ref){
                    list[i]=i1ref;
                }
                if(list[i]==i2ref) list[i]=list[i]-1;
            }        //70 continue

            howmany[i1ref]=howmany[i1ref]+howmany[i2ref];

            for (k=(i2ref+1);k<=kk;k++){
                howmany[k-1]=howmany[k];
            }

          }    //end for each k
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

        if(intParam==1){
            for (k=k1;k>=k2;k--){
                if (CHr[k]>=CHr_max){
                    CHr_group=k;
                    CHr_max=CHr[k];

                    //Pour évider les clusters vides
                    realk=0;
                    unique=0;
                    CHk = 0;

                    for(int i=1; i<=k; i++){
                        unique=0;
                        for(int j=1; j<=n; j++){
                            if(listr[CHr_group][j]==i && unique==0){
                                CHk++;
                                CH_conversion[i] = CHk;
                                unique=1;
                                realk++;
                            }
                        }
                    }

                    for(int iz=1; iz<=N; iz++){
                        listr[CHr_group][iz]=CH_conversion[listr[CHr_group][iz]];
                        Strouve[iz-1]=CH_conversion[listr[CHr_group][iz]];
                    }
                    CHr_group=realk;
                }
            }
        }else if(intParam==2){
            for (k=k1;k>=k2;k--){
                if (Wr[k]<=W_min){
                    //pour connaitre le nombre de partition adéquate.
                    W_group=k;
                    W_min=Wr[k];

                    realk=0;
                    unique=0;
                    wk = 0;

                    for(int i=1; i<=k; i++){
                        unique=0;
                        for(int j=1; j<=n; j++){
                            if(listr[W_group][j]==i && unique==0){
                                wk++;
                                W_conversion[i] = wk;
                                unique=1;
                                realk++;
                            }
                        }
                    }

                    for(int iz=1; iz<=N; iz++){
                        listr[W_group][iz]=W_conversion[listr[W_group][iz]];
                        Strouve[iz-1]=W_conversion[listr[W_group][iz]];
                    }
                    W_group=realk;
                }

            }
        }


    }  //fin random start

    // Print results
    switch (intParam){
        case 1:
        {
            strcpy(criteria, "CH");
            conv2sameRef(Strouve,Sref,N);
            outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr,CHr,k1,k2,monTableau);
        }break;

        case 2:
        {
            strcpy(criteria, "BH");
            conv2sameRef(Strouve,Sref,N);
            outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,W_group,W_max,listr,Wr,k1,k2,monTableau);
        }break;

    }

    // End timer
    tend2=time(NULL);                // get the current calendar time

    // Compute execution time
    texec2=difftime(tend2,tbegin2);    // tend-tbegin (result in second)
    fprintf (Output4,"%.3f;\n",texec2);

    //}// fin de random start

    // Print results


    // *********************Close output files ***************************
    fclose(Output4);
    //*********************** Remove matrix ******************************

    for (i=0;i<=kmax;i++){
        delete [] sx[i];
        delete [] sx2[i];
        delete [] xbar[i];
        delete [] var[i];
        delete [] listr[i];
        delete [] howmanyr[i];
    }

    delete [] sx;
    delete [] sx2;
    delete [] xbar;
    delete [] var;
    delete [] listr;
    delete [] howmanyr;

    delete [] Dvec;

    delete [] CHr;
    delete [] Wr;
    delete [] Wr_ln;
    delete [] diff_W;
    delete [] V_W;

    delete [] SSEr;

    delete [] vect;
    delete [] mean;
    delete [] weight;

    delete [] list;
    delete [] no;
    delete [] iordre;

    delete [] howmany;
    delete [] nobest;
    delete [] nnitr;

    delete [] ishort;

    delete [] nameb;
    delete [] nk;
    delete [] distances_RF_norm;

    for (int i=0;i<n;i++)
    {
        delete [] tree_cluster_leaves[i];
    }
    delete [] tree_cluster_leaves;


    return 0;
}


 //      end
//************************End of Main

//******************************************************************************
//**********************************FUNCTIONS***********************************
//******************************************************************************

void ReadData1(int &n,int &nmax,int &p,int &pmax,double** mat/* ,double* coord */,int* ishort,double* weight,double* colsum,int &ntran, char* nameb, int N)
{
    int p1=0,p2=0;

    int j=0;
    int   nmat=2; // (orientation of data))

    //Read matrix parameters
    n = N;
    p = N;
    //printf("\nData:\nn:%d p:%d\n", n,p);

    if(n>nmax)    //if(n.gt.nmax) Stop
    {
        printf ("Too many objects. Use a sample of objects or recompile program to increase nmax.");                //     +'Too many objects. Use a sample of objects or recompile program.'
        exit(1);
    }

    if(p>pmax)    // if(p.gt.pmax) Stop 'Too many variables. Recompile the program.'
    {
        printf ("Too many variables. Use a sample of objects or recompile program to increase pmax.");                //     +'Too many objects. Use a sample of objects or recompile program.'
        exit(1);
    }


   switch (nmat)   //goto (10,14,18) nmat
   {
       case 1:
       {
            break;//goto 22;
       }
       case 2:
       {
            break;//goto 22;
       }

    // To read the file of QTC variables:
       case 3:
       {
            printf ("How many position (p1) and QTC variables (p2)?\n");
            printf ("E.g., p1 = 0 or 1 or 2; p2 = 166.  Answer with 2 numbers:");
            scanf("%d %d",&p1, &p2);

            if(p2>pmax){
                printf ("Too many variables. Recompile the program to increase pmax.");
                exit(1);
            }

            p=p2;
            printf ("\n");

            break;
        }
    }


   //fclose(Input1);

    for (j=1;j<=p;j++){
        ishort[j]=j;
        weight[j]=1.0;
    }

    strcpy(nameb,"../output/stat.csv");
    if((Output4 = fopen(nameb,"a"))==NULL){
        printf("\n%s: result file open failed...",nameb);
        exit(1);
    }

}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

void Assign(int &iran,int &n,int &nmax,int &k1,int &kmax,int* list,int* howmany,int* no,int &idebug,int &iassign,int &iseed, int random_number)

{
    int k=0, i=0, ii=0, kk=0, how=0, isum=0;
    char namea[MAX_PATH_LENGTH];
    double turn=0;

    if ((iassign==1) || (iassign==2)){
        how=n/(k1*1.0);
        for (k=1;k<=(k1-1);k++) {howmany[k]=how;}
        howmany[k1]=n-(k1-1)*how;
        ii=0;

        for (k=1;k<=k1;k++){
            for (kk=1;kk<=howmany[k];kk++){
               ii++;
               list[ii]=k;
            }
        }


        if(iassign==1) return;
        // Assign objects at random to the groups
        if(iran==1){
            for (i=1;i<=(random_number+100);i++)  turn=rand()/(1.0*(rand() % RAND_MAX_VALUE));
        }                            //end if
        Permute(iseed,n,nmax,list);
        return;
    }else if (iassign==3){
        // Read file of group assignments.
        // First line: how many objects in each group?
        // Then, read members of each group on one line (list of object numbers).

        printf ("Name of file of group assignments?");        //60 write(*,*) 'Name of file of group assignments?'
        scanf ("%s",namea);        //read(*,*) namea

        FILE *Input3;
        if ((Input3 = fopen(namea,"r"))==0) { printf("\n %s :Open Failed....",namea); exit(1); }
        printf ("File of group assignments: %s\n",namea);

        for (k=1;k<=k1;k++){
            fscanf(Input3,"%d",&howmany[k]);
        }

        isum=0;
        for (k=1;k<=k1;k++){
            isum=isum+howmany[k];
        }

        if(isum!=n){
            printf("Objects assigned to groups do not sum to n.");
            exit(1);
        }

        for (i=1;i<=n;i++) {
            list[i]=-1;
        }

        for (k=1;k<=k1;k++){
            for (i=1;i<=howmany[k];i++){
                fscanf(Input3, "%d", &no[i]);
            }
            for (i=1;i<=howmany[k];i++){
                list[no[i]]=k;
            }
        }

        for (i=1;i<=n;i++){
            if(list[i]==-1){
                printf("Overlapping assignments to groups.");
                exit(1);
            }
        }
        fclose(Input3);
        return;
    }else{
        printf("Wrong perameter <iassign> in function <Assign>.");
        exit(1);
    }

}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

void CompSST(int &n,int &nmax,int &p,int &pmax,double** mat,double* weight,int* ishort,double &SST)
{
    double    sx=0,sx2=0,var=0,temp=0,dfln=0;     //Real*8 mat(nmax,pmax),weight(pmax),sx,sx2,var,temp,dfln,SST
    int j=0, i=0;
    dfln=n;        //dfln=dfloat(n)
    SST=0.0;                //SST=0.0

    for (j=1;j<=p;j++)        // do 22 j=1,p
    {
        sx=0.0;        //sx=0.0
        sx2=0.0;        //sx2=0.0
        for (i=1;i<=p;i++)        // do 20 i=1,n
        {
            temp=mat[i-1][j-1];
            sx=sx+temp;                //sx=sx+temp
            sx2=sx2+temp*temp;        //20 sx2=sx2+temp*temp
        }
        var=sx2-(sx*sx/(1.0*dfln));            //var=sx2-(sx*sx/dfln)
        SST=SST+var*weight[ishort[j]];        //22 SST=SST+var*weight(ishort(j))
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

void Permute(int &iseed,int &n,int &nmax,int *iordre)
{
    int i=0,j=0,itemp=0,km1=0,m=0;            //Integer iseed,n,i,j,itemp,km1,m
    //Integer iordre(nmax)
    m=n;        //m=n
    km1=n-1;    //km1=n-1

    for (i=1;i<=km1;i++)            //      do 10 i=1,km1
    {
m8:        // j = 1 + rand(iseed)*m;        //    8    j = 1 + rand(iseed)*m
        j=1+(rand()/RAND_MAX_VALUE)*m;

        if(j>m) goto m8;            //if(j.gt.m) goto 8
        itemp = iordre[m];            //itemp = iordre(m)
        iordre[m] = iordre[j];        //iordre(m) = iordre(j)
        iordre[j] = itemp;            //iordre(j) = itemp
        m = m-1;                    //m = m-1
    }            //   10 continue
    return;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

//compute rand index
double f_RI(int Strouve[],int Sref[],int N)
{
    double comb = 1.0;

    for (int i=N; i>=(N-2+1); i--)
    {
        comb*=i;
    }

    comb/=2.0;

    double a=0.0;
    double b=0.0;

    for (int i=0; i<N-1; i++){
        for (int j=i+1; j<N; j++){
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
double f_ARI(int Strouve[],int Sref[],const char *K_real,int group,int N)
{
    int kReal = atoi(K_real);
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


    for(int i=0;i<N;i++){
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
    double d=0.0;

    double comb = 1.0;

    for (int i=N; i>=(N-2+1); i--)
    {
        comb*=i;
    }

    comb/=2.0;

    for (int i=0; i<N-1; i++){
        for (int j=i+1; j<N; j++){
            if(Sref[i]!=Sref[j]){
                if(Strouve[i]!=Strouve[j]){
                    d++;
                }else{
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
void outStat(int Strouve[],int Sref[],char *criteria,int N,char *N_especes,char *percent,const char *K_real,int group/* ,double RI,double ARI */,double score,int **listr,double *allScore,int k1, int k2, vector <string> monTableau)
{

    //Compute Rand index between Strouve and Sref
    double RI = f_RI(Strouve,Sref,N);

    //Compute Rand index Adjusted between Strouve and Sref
    double ARI = f_ARI(Strouve,Sref,K_real,group,N);

    fprintf (Output4,"%s;",criteria);
    fprintf (Output4,"%i;",N);
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
    for (int p=1; p<=N; p++){
        if(p==N){
            //fprintf (Output4,"%i%s",listr[group][p]," ");
            fprintf (Output4,"%i%s",Strouve[p-1]," ");
        }else{
            //fprintf (Output4,"%i%s",listr[group][p]," <> ");
            fprintf (Output4,"%i%s",Strouve[p-1]," <> ");
        }

    }
    fprintf (Output4,");");

    
    // for(int i=k2; i<=k1; i++){
        // cout<<"K = "<<i<<" "<<criteria<<" = "<<allScore[i]<<" : ";
        // for(int j=1; j<=N; j++){
            // cout<<listr[i][j]<<" <> ";
        // }
        // cout<<endl;
    // }
    
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
        for (int p=1; p<=N; p++){
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
    
    cout<<"1) stat.csv - for clustering statistics;"<<endl;
    cout<<"2) output.txt - for cluster content."<<endl;
}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double FO_super_tree(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, int alpha)
{
    double *clusterK_same = new double [kmax+1];
    int *nk_CH = new int [kmax+1];
    int cluster_k=0;
    double RF = 0.0;
    double Dref=0;     //Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
    int    kref=0;        //Integer list(nmax),howmany(kmax),kref
    //Integer ishort(pmax)
    // Compute squared distances to group centroids. Assign objects to nearest one
    SSE=0;        //SSE=0.0

    int k_source = 0;
    int new_k = 0;
    int old_k = 0;
    int nb_cluster_dest = 0;
    int nb_cluster_source = 0;
    double FO_old = 0.0;
    double FO_new = 0.0;
    double tmp_calc_dest = 0.0;
    double tmp_calc_source = 0.0;

    //initialisation des variables
    for(int k=1;k<=kmax; k++){
        nk_CH[k]=0;
        clusterK_same[k]=0.0;
    }

    //Compter le nombre d'éléments par cluster
    for(int k=1;k<=kmax; k++){
        nk_CH[list[k]]++;
    }
    
    //Créer l'arbre consensus pour chaque cluster
    if(withConsensus == true){
        vector <string> listeConsensusCluster;
        double *distances = new double[6];
        for (int j=0; j<4; j++){
            distances[j]=0.0;
        }
        distances[5] = alpha;
        for(int k=1;k<=kmax; k++){  //pour chaque cluster
            int indice_arbre= 0; //indice de l'arbre
            //créer le consensus
            string consensek = "init";
            
            //ajouter le consensus à listeConsensusCluster
            listeConsensusCluster[k - 1] =  consensek;
            for (int i = 0; i < nk_CH[list[k]]; i++){   //pour chaque arbre du cluster k
                //Calcul des distances entre chaque abre du cluster et le consensus
                main_hgt(consensek, monTableau[indice_arbre], distances);
                clusterK_same[k] += distances[0];
            }
        }
    }
    else{
        //compute for each cluster initially, SSW value (intra groupe distance)
        //compute SSW
        for (int i=1;i<n;i++){
            cluster_k=list[i];
            for (int j=i+1;j<=n;j++){
                if (list[j]==cluster_k){
                    
                    RF = mat[i-1][j-1];
                    clusterK_same[cluster_k]+=RF;
                }
            }
        }
    }

    for (int k=1;k<=kk;k++){
        if(nk_CH[k]>1){
            FO_old += ((clusterK_same[k])/(1.0*nk_CH[k]));
        }
    }
    for (int i=1;i<=n; i++){
        if(nk_CH[list[i]]>1){
            for (int k=1;k<=kk;k++)            //do 12 k=1,kk
            {
                //Calcul de la distance RF de chaque point i
                // et assignation du point i au bon cluster
                // Compute a RF distance to the centroid k

                //test si le point i n'appartenait pas initiallement à k
                //k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
                if((list[i]!=k) && (nk_CH[list[i]]>1)){
                    //Pour le cluster source
                    k_source = list[i];
                    nb_cluster_source = nk_CH[k_source];

                    tmp_calc_source = clusterK_same[k_source];
                    nb_cluster_source -=1;
                    if (nb_cluster_source==1){
                        tmp_calc_source = 0.0;
                    }else{
                        for(int j=1;j<=n; j++){
                            if(list[j]==k_source && i!=j){
                                tmp_calc_source -= mat[i-1][j-1];
                            }
                        }
                    }

                    nb_cluster_dest = nk_CH[k];

                    if((nk_CH[k]>1) && (nk_CH[list[i]]>1)){
                        FO_new = FO_old - (clusterK_same[k]/(1.0*nk_CH[k]));
                        FO_new = FO_new - (clusterK_same[list[i]]/(1.0*nk_CH[list[i]]));
                    }
                    else if(nk_CH[list[i]]>1){
                        FO_new = FO_old - (clusterK_same[list[i]]/(1.0*nk_CH[list[i]]));
                    }
                    else if(nk_CH[k]>1){
                        FO_new = FO_old - (clusterK_same[k]/(1.0*nk_CH[k]));
                    }
                    else{
                        FO_new = FO_old;
                    }

                    //compute Function objective
                    //Pour le cluster de destination
                    tmp_calc_dest = clusterK_same[k];
                    for(int j=1;j<=n; j++){
                        if(list[j]==k){
                            tmp_calc_dest += mat[i-1][j-1];
                        }
                    }
                    nb_cluster_dest +=1;
                    if(nb_cluster_dest>=1){
                        FO_new = FO_new + (tmp_calc_dest/(1.0*nb_cluster_dest));
                    }

                    //Pour le cluster source
                    if(nb_cluster_source>=1){
                        FO_new = FO_new + (tmp_calc_source/(1.0*nb_cluster_source));
                    }

                    if(FO_new<FO_old){
                        Dref=FO_new;
                        kref=k;

                        //mise à jour de nk_CH[]
                        nk_CH[k] = nb_cluster_dest;
                        nk_CH[list[i]] = nb_cluster_source;

                        //mise à jour de la distance intra-groupe des deux clusters modifiés
                        clusterK_same[list[i]] = tmp_calc_source;
                        clusterK_same[k] = tmp_calc_dest;

                        howmany[k] = nb_cluster_dest;
                        howmany[list[i]] = nb_cluster_source;

                        SSE=SSE+Dref;         //SSE=SSE+Dref

                        //mise à jour de la fonction objective FO_old
                        FO_old = FO_new;

                        //mise à jour la liste de distribution des elements
                        list[i] = k;

                        //A VOIR SI UTILE
                        new_k = k;
                        old_k = list[i];

                    }else{
                        Dref=FO_old;
                    }
                }
            }
        }
    }

    delete [] clusterK_same;
    delete [] nk_CH;

  return Dref;

}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double DistanceCH(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new,double facteur) {

    double SSB=0.0;
    double SSW=0.0;
    double dist_all = 0.0;
    double RF;
    double distance_total = 0.0;
    int *nk_CH = new int [kmax+1];
    int k_cluster = 0;

    for(int k=1;k<=kmax; k++){
        nk_CH[k]=0;
    }

    for(int k=1;k<=kmax; k++){
        nk_CH[list[k]]++;
    }

    for(int k=1;k<=kmax; k++){
        if(nk_CH[k]!=0){
            k_cluster++;
        }
    }

    //compute dist_all
    for (int i=1;i<n;i++){
        for (int j=i+1;j<=n;j++){
            RF = mat[i-1][j-1];
            dist_all += RF;
        }
    }

    dist_all = dist_all/(1.0*n);

    //compute SSW
    SSW = FO_new;

    //compute SSB
    SSB = dist_all - SSW;

     if((fabs(SSW)>0.000001) && (k_cluster>1)){
        distance_total=(SSB/(1.0*SSW))*((1.0*n-k_cluster)/((1.0*k_cluster)-1.0));
    }

    if(fabs(SSW)<=0.000001  && (k_cluster>1)){
        distance_total=10000000.0*SSB*((n-k_cluster)/((1.0*k_cluster)-1.0));
    }

    delete [] nk_CH;

    return distance_total;

}

// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

double FO_W(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau)
{
    double *clusterK_same = new double [kmax+1];
    int *nk_W = new int [kmax+1];
    int cluster_k=0;
    double RF = 0.0;
    double Dref=0;       //Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
    int    kref=0;        //Integer list(nmax),howmany(kmax),kref
    //Integer ishort(pmax)
    // Compute squared distances to group centroids. Assign objects to nearest one
    SSE=0;        //SSE=0.0

    int k_source = 0;
    int new_k = 0;
    int old_k = 0;
    int nb_cluster_dest = 0;
    int nb_cluster_source = 0;
    double FO_old = 0.0;
    double FO_new = 0.0;
    double tmp_calc_dest = 0.0;
    double tmp_calc_source = 0.0;

    //initialisation des variables
    for(int k=1;k<=kmax; k++){
        nk_W[k]=0;
        clusterK_same[k]=0.0;
    }

    //Compter le nombre d'éléments par cluster
    for(int k=1;k<=kmax; k++){
        nk_W[list[k]]++;
    }

    //compute for each cluster initially, SSW value (intra groupe distance)
    //compute SSW
    for (int i=1;i<n;i++){
        cluster_k=list[i];
        for (int j=i+1;j<=n;j++){
            if (list[j]==cluster_k){
                RF = mat[i-1][j-1];
                clusterK_same[cluster_k]+=RF;
            }
        }
    }

    for (int k=1;k<=kk;k++){
        if(nk_W[k]>1){
            FO_old += ((2.0*clusterK_same[k])/(1.0*nk_W[k]*(nk_W[k]-1)));
        }
    }

    if(kk==1){
        Dref=FO_old;
    }else{
        for (int i=1;i<=n; i++){
            if(nk_W[list[i]]>1){
                for (int k=1;k<=kk;k++){
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
                            for(int j=1;j<=n; j++){
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
                        for(int j=1;j<=n; j++){
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

                            SSE=SSE+Dref;         //SSE=SSE+Dref

                            //mise à jour de la fonction objective FO_old
                            FO_old = FO_new;

                            //mise à jour la liste de distribution des elements
                            list[i] = k;

                            //A VOIR SI UTILE
                            new_k = k;
                            old_k = list[i];
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

double DistanceW(int &n,int &kmax,double** mat, int* list, double** Ww, double FO_new, double facteur) {


    double distance_total = 100000000.0;
    double *clusterK_same = new double [kmax+1];
    int *nk_W = new int [kmax+1];
    int k_cluster = 0;


    for(int k=1;k<=kmax; k++){
        nk_W[k]=0;
        clusterK_same[k]=0.0;
    }


    for(int k=1;k<=kmax; k++){
        nk_W[list[k]]++;
    }

    for(int k=1;k<=kmax; k++){
        if(nk_W[k]!=0){
            k_cluster++;
        }
    }

     if(k_cluster!=n){
        // distance_total=(FO_new/(1.0*(n-k_cluster)));
        distance_total=(FO_new/k_cluster);
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
void conv2sameRef(int *Strouve,int *Sref, int n)
{
    int k = 0;

    cout<<"Number of trees in the input file: "<<n<<endl;
    cout<<"Partition found: "<<endl;
    for(int i=0; i<n; i++){
        cout<<Strouve[i]<<" <> ";
    }
    cout<<endl;

    // cout<<"Sref"<<endl;
    // for(int i=0; i<n; i++){
        // cout<<Sref[i]<<" <> ";
    // }
    // cout<<endl;

    //To know the number of cluster
    for(int i=0; i<n; i++){
        if(Strouve[i]>k){
            k=Strouve[i];
        }
    }

    cout<<"Number of clusters (K) found: "<<k<<endl;

    int *nk_trouve = new int [k];
    int *nk_ref = new int [k];

    for(int i=0;i<k; i++){
        nk_trouve[i]=0;
        nk_ref[i]=0;
    }


    //Compter le nombre d'éléments par cluster
    //For Strouve
    for(int i=0;i<k; i++){
        nk_trouve[Strouve[i]]++;
    }

    //For Sref
    for(int i=0;i<k; i++){
        nk_ref[Sref[i]]++;
    }

    delete []nk_trouve;
    delete []nk_ref;
}
