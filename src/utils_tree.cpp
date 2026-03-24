

//  utils_tree.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

/* =========================================================================
 * utils_tree.cpp — Fonctions utilitaires pour le traitement d'arbres
 * phylogénétiques dans le projet K-Means Phylogenetic Trees Clustering.
 *
 * Ce fichier regroupe :
 *   - Calcul de distances entre arbres (Robinson-Foulds, bipartition, LS)
 *   - Algorithme de Floyd-Warshall pour les plus courtes distances
 *   - Construction et parcours de structures d'arbres (TNoeud, DescTree)
 *   - Lecture et écriture au format Newick et matriciel
 * ========================================================================= */

#include "utils_tree.hpp"

/**
 * Calcule la distance de bipartition entre deux ensembles de bipartitions
 * B et B1 pour un arbre de n feuilles.
 *
 * Pour chaque bipartition non triviale (ni singleton ni complémentaire)
 * de B, on cherche la bipartition la plus proche dans B1 (et réciproquement).
 * Le résultat est la moyenne des sommes de discordances minimales dans les
 * deux sens.
 *
 * Paramètres :
 *   B, B1 - tableaux de bipartitions des deux arbres
 *   n     - nombre de feuilles
 *
 * Retourne la distance de bipartition moyenne.
 */
double BipartitionDistance (int **B, int ** B1,int n)
{
    int i,j,k,k1,cpt1,cpt2,cpt3,cpt4;
    double c1,c2,min1 = 0.0,min2 = 0.0,cptB=0.0,cptB1=0.0;
    int *flag,*flag1,**Bi,**Bi1;
    Bi = safe_malloc<int*>(2*n-2);
    Bi1 = safe_malloc<int*>(2*n-2);
    for (i=0;i<2*n-2;i++)
    {
        Bi[i] = safe_malloc<int>(2*n);
        Bi1[i] = safe_malloc<int>(2*n);
        if ((Bi1[i]==NULL)||(Bi[i]==NULL))
        {
            printf(" Data matrix is too large\n ");
            exit(2);
        }
    }
    flag = safe_malloc<int>(2*n);
    flag1 = safe_malloc<int>(2*n);

    /* Marquer les bipartitions triviales (singleton ou complémentaire) */
    for(i=1;i<=2*n-3;i++){
        flag[i] = flag1[i] = 0;
        for(j=1;j<=n;j++){
            flag[i] = flag[i] + B[i][j];
            flag1[i] = flag1[i] + B1[i][j];
        }
        if(flag[i] == 1 || flag[i] == n-1)
            flag[i] = 1;
        else
            flag[i] = 0;
        if(flag1[i] == 1 || flag1[i] == n-1)
            flag1[i] = 1;
        else
            flag1[i] = 0;
    }

    /* Copier uniquement les bipartitions non triviales dans Bi et Bi1 */
    k=k1=1;
    for(i=1;i<=2*n-3;i++){
        if(flag[i] == 0){
            for(j=1;j<=n;j++)
                Bi[k][j] = B[i][j];
            k++;
        }
        if(flag1[i] == 0){
            for(j=1;j<=n;j++)
                Bi1[k1][j] = B1[i][j];
            k1++;
        }
    }

    /* Calculer la discordance minimale dans les deux sens */
    for(i=1;i<=n-3;i++){
        for(j=1;j<=n-3;j++){
            cpt1=cpt2=cpt3=cpt4=n;
            for(k=1;k<=n;k++){
                if(Bi[i][k]==Bi1[j][k])
                    cpt1--;
                else
                    cpt2--;
                if(Bi1[i][k]==Bi[j][k])
                        cpt3--;
                    else
                        cpt4--;
            }
            c1 = ((double)cpt1<(double)cpt2)?(double)cpt1:(double)cpt2;
            if(j==1 || c1 < min1)
                min1 = c1;
            c2 = ((double)cpt3<(double)cpt4)?(double)cpt3:(double)cpt4;
            if(j==1 || c2 < min2)
                min2 = c2;
        }
        cptB  += min1;
        cptB1 += min2;
    }

    for (i=0;i<2*n-2;i++)
    {
        free(Bi[i]); free(Bi1[i]);
    }

    free(Bi);
    free(Bi1);
    free(flag); free(flag1);

    return (cptB + cptB1) / 2.0;
}


/**
 * Exécute l'algorithme de Floyd-Warshall sur la matrice d'adjacence
 * pour calculer toutes les distances de plus court chemin.
 *
 * Le résultat est écrit dans DIST. Le paramètre kt représente le nombre
 * de branches terminales à ignorer (les n premiers sommets sont les feuilles).
 *
 * Paramètres :
 *   Adjacence - matrice d'adjacence de l'arbre
 *   DIST      - matrice de distances résultante
 *   n         - nombre d'espèces (feuilles)
 *   kt        - nombre de branches terminales ignorées
 */
void Floyd(double ** Adjacence , double ** DIST,int n,int kt)
{
    int i,j,k;

    for(i=1;i<=2*n-2-kt;i++) {
        for(j=1;j<=2*n-2-kt;j++) {
            if(i==j) {
                DIST[i][j] = 0;
            } else {
                DIST[i][j] = Adjacence[i][j];
            }
        }
    }

    for(i=1;i<=2*n-2-kt;i++) {
        for(j=1;j<=2*n-2-kt;j++) {
            for(k=1;k<=2*n-2-kt;k++) {
                if((DIST[j][i] + DIST[i][k]) < DIST[j][k]) {
                    DIST[j][k] = DIST[j][i] + DIST[i][k];
                }
            }
        }
    }
}


/**
 * Variante de Floyd-Warshall remplissant deux matrices de distances
 * simultanément (DIST et DIST2) à partir de la même matrice d'adjacence.
 *
 * Permet de conserver une copie initiale des distances sans écraser la
 * première. Voir Floyd() ci-dessus pour la description des paramètres.
 */
void Floyd(double ** Adjacence , double ** DIST,double **DIST2,int n,int kt)
{
    int i,j,k1;

    for(i=1;i<=2*n-2-kt;i++)
        for(j=1;j<=2*n-2-kt;j++)
        {
            if(i==j)
                DIST[i][j] = DIST2[i][j] = 0;
            else
                DIST[i][j] = DIST2[i][j] = Adjacence[i][j];
        }

    for(i=1;i<=2*n-2-kt;i++) {
        for(j=1;j<=2*n-2-kt;j++) {
            for(k1=1;k1<=2*n-2-kt;k1++) {
                if((DIST[j][i] + DIST[i][k1]) < DIST[j][k1]) {
                    DIST[j][k1] = DIST2[j][k1] = DIST[j][i] + DIST[i][k1];
                }
            }
        }
    }
}


/**
 * Reconstruit la matrice d'adjacence à partir des listes d'arêtes (ARETE)
 * et de leurs longueurs (LONGUEUR).
 *
 * Toutes les entrées sont d'abord initialisées à INFINI, puis les arêtes
 * connues sont renseignées symétriquement. Les branches terminales (kt)
 * sont exclues de la reconstruction.
 *
 * Paramètres :
 *   Adjacence - matrice à remplir (modifiée en place)
 *   ARETE     - tableau des paires de sommets formant les arêtes
 *   LONGUEUR  - longueurs associées à chaque arête
 *   size      - nombre d'espèces
 *   kt        - nombre de branches terminales ignorées
 */
void loadAdjacenceMatrix( double **Adjacence, long int *ARETE, double *LONGUEUR,int size,int kt){

    int i,j;

    /* Initialisation de toutes les distances à l'infini */
    for(i=1;i<=2*size-2;i++)
        for(j=1;j<=2*size-2;j++){
            Adjacence[i][j] = Adjacence[j][i] = INFINI;
        }

    /* Remplissage symétrique des arêtes connues */
    for(i=1;i<=2*size-3-kt;i++){
        Adjacence[ARETE[2*i-2]][ARETE[2*i-1]] = LONGUEUR[i-1];
        Adjacence[ARETE[2*i-1]][ARETE[2*i-2]] = LONGUEUR[i-1];
    }
}


/**
 * Calcule un ordre circulaire approximatif (outward permutation) sur
 * la matrice de distances D.
 *
 * L'ordre est stocké dans X[1..n]. Les extrémités i1 et j1 sont fixées
 * comme premier et dernier élément. L'algorithme insère les éléments
 * restants un par un en minimisant le critère de circularité.
 *
 * Paramètres :
 *   D       - matrice de distances
 *   X       - tableau résultat de l'ordre circulaire
 *   i1, j1  - indices des extrémités fixées
 *   n       - nombre d'éléments
 */
void odp1(double **D, int *X, int *i1, int *j1, int n)
{
    double S1,S;
    int i,j,k = 0,a,*Y1;

    Y1 = safe_malloc<int>(n+1);

    for(i=1;i<=n;i++)
        Y1[i]=1;

    X[1]=*i1;
    X[n]=*j1;
    if (n==2){
        free(Y1);
        return;
    }
    Y1[*i1]=0;
    Y1[*j1]=0;
    for(i=0;i<=n-3;i++)
    {
        a=2;
        S=0;
        for(j=1;j<=n;j++)
        {
            if (Y1[j]>0)
            {
                S1= D[X[n-i]][X[1]]-D[j][X[1]]+D[X[n-i]][j];
                if ((a==2)||(S1<=S))
                {
                    S=S1;
                    a=1;
                    X[n-i-1]=j;
                    k=j;
                }
            }
        }
        Y1[k]=0;
    }
    free(Y1);
}


/**
 * Détermine les arêtes d'un arbre à partir de la matrice de distances DI.
 *
 * Les arêtes sont stockées dans ARETE (paires de sommets) et leurs longueurs
 * dans LONGUEUR. L'algorithme construit l'arbre en insérant les feuilles
 * dans l'ordre circulaire fourni par odp1, en décomposant chaque chemin.
 *
 * Si binaire == 0, les branches internes quasi-nulles (≤ 2*epsilon_ba) sont
 * contractées de manière itérative pour produire un arbre non binaire.
 *
 * Paramètres :
 *   DI      - matrice de distances des feuilles
 *   ARETE   - tableau des arêtes résultant (modifié)
 *   LONGUEUR - longueurs des arêtes (modifié)
 *   n       - nombre de feuilles
 *   binaire - 1 = conserver l'arbre binaire, 0 = contracter les branches nulles
 *
 * Retourne le nombre de branches terminales contractées (kt).
 */
int Tree_edges (double **DI, long int *ARETE, double *LONGUEUR, int n,int binaire)
{

    struct EDGE { unsigned int U; unsigned int V; double LN;};
    struct EDGE *Path,*Tree;
    int i,j,k,p,P,*X;
    double S,DIS,DIS1,*L,**D;
    int pasfini = 1;
    int kt=0;
    long SomToDel,OtherSom;

    X = safe_malloc<int>(n+1);
    L = safe_malloc<double>(n+1);
    Tree = safe_malloc<struct EDGE>(2*n-2);
    Path = safe_malloc<struct EDGE>(n+2);

    D=(double **) malloc((n+1)*sizeof(double*));

    for (i=0;i<=n;i++)
    {
        D[i]=(double*)malloc((n+1)*sizeof(double));

        if (D[i]==NULL)
        {
            printf("Data matrix is too large"); exit(1);
        }
    }

    i=1; j=n;
    odp1(DI,X,&i,&j,n);

    for (i=1;i<=n;i++)
    {
        for (j=1;j<=n;j++)
            D[i][j]=DI[i][j];
    }

    /* Vérification de l'équivalence des topologies */
    L[1]=D[X[1]][X[2]];
    Path[1].U=X[1];
    Path[1].V=X[2];

    p=0;
    P=1;

    for(k=2;k<=n-1;k++)
    {

        DIS=(D[X[1]][X[k]]+D[X[1]][X[k+1]]-D[X[k]][X[k+1]])/2;
        DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;

        S=0.0;
        i=0;

        if (DIS>2*epsilon_ba)
        {
            while (S<DIS-(epsilon_ba))
            {
                i=i+1;
                S=S+L[i];
            }
        }
        else { DIS=0; i=1; }

        Tree[p+1].U=n+k-1;
        Tree[p+1].V=Path[i].V;
        Tree[p+1].LN=S-DIS;
        if (Tree[p+1].LN<epsilon_ba) Tree[p+1].LN=2*epsilon_ba;

        for (j=i+1;j<=P;j++)
        {
            Tree[p+j-i+1].U=Path[j].U;
            Tree[p+j-i+1].V=Path[j].V;
            Tree[p+j-i+1].LN=L[j];
            if (L[j]<2*epsilon_ba) L[j]=2*epsilon_ba;
        }
        p=p+P-i+1;

        Path[i].V=n+k-1;
        Path[i+1].U=n+k-1;
        Path[i+1].V=X[k+1];
        L[i]=L[i]+DIS-S;
        L[i+1]=DIS1;
        P=i+1;
    }

    for (i=1;i<=P;i++)
    {
        Tree[p+i].U=Path[i].U;
        Tree[p+i].V=Path[i].V;
        Tree[p+i].LN=L[i];
    }

    for (i=1;i<=2*n-3;i++)
    {
        if (fabs(Tree[i].LN-epsilon_ba)<=2*epsilon_ba)
            Tree[i].LN=0.0;
        ARETE[2*i-2]=Tree[i].U;
        ARETE[2*i-1]=Tree[i].V;
        LONGUEUR[i-1]=Tree[i].LN;
        if (LONGUEUR[i-1]<2*epsilon_ba) LONGUEUR[i-1] = 2*epsilon_ba;
    }

    /* Contraction des branches internes quasi-nulles (ajout 22 avril 2005) */
    while((pasfini==1)&&(binaire==0)){
        pasfini = 0;
        for (i=1;i<=2*n-3-kt;i++){
            if((LONGUEUR[i-1] == 2*epsilon_ba)&&(ARETE[2*i-2]>n)&&(ARETE[2*i-1]>n)){
                /* branche interne de longueur quasi-nulle : contraction */
                if(ARETE[2*i-2] > ARETE[2*i-1]){
                    SomToDel = ARETE[2*i-2];
                    OtherSom = ARETE[2*i-1];
                }
                else{
                    SomToDel = ARETE[2*i-1];
                    OtherSom = ARETE[2*i-2];
                }
                pasfini=1;

                for (j=1;j<=2*n-3-kt;j++){

                    if(j!=i){
                        if(ARETE[2*j-2] == SomToDel){
                            ARETE[2*j-2] = OtherSom;
                        }
                        else if(ARETE[2*j-1] == SomToDel){
                            ARETE[2*j-1] = OtherSom;
                        }
                    }
                }
                for(j=i;j<=2*n-3-kt;j++){
                    LONGUEUR[j-1] = LONGUEUR[j];
                    ARETE[2*j-1] = ARETE[2*(j+1)-1];
                    ARETE[2*j-2] = ARETE[2*(j+1)-2];
                }
                for(j=1;j<=2*n-3-kt-1;j++){
                    if(ARETE[2*j-2] > SomToDel)
                        ARETE[2*j-2]--;
                    if(ARETE[2*j-1] > SomToDel)
                        ARETE[2*j-1]--;
                }
                break;
            }
        }
        if(pasfini)
            kt++;
    }
    free(X);
    free(Tree);
    free(L);
    free(Path);

    for (i=0;i<=n;i++)
        free(D[i]);

    free(D);

    return kt;
}


/**
 * Construit la table de bipartitions d'un arbre à partir de sa matrice de
 * distances D, en utilisant l'ordre circulaire de odp1.
 *
 * Chaque ligne de B représente une bipartition de l'arbre : B[i][j] vaut 1
 * si l'espèce j est dans le sous-ensemble gauche de la bipartition i.
 * PLACE[i] contient l'indice de la i-ème bipartition dans l'ordre de
 * la colonne maximale.
 *
 * Paramètres :
 *   D     - matrice de distances des feuilles
 *   B     - tableau de bipartitions résultant (modifié)
 *   PLACE - tableau d'indices de tri (modifié)
 *   n     - nombre de feuilles
 *
 * Retourne le nombre total de bipartitions non triviales (m).
 */
int Bipartition_Table (double **D, int **B, int *PLACE, int n)
{
    int i,j,k,l,l1,*MaxCol,*X,EdgeNumberPath,m,uv = 0,PlaceNumber,edge,*Path,M,F;
    double S,DIS,DIS1,*LengthPath;
    double EPS=1.e-5;
    double EPS1=1.e-2;

    MaxCol=(int *)malloc((2*n)*sizeof(int));
    X=(int *)malloc((2*n+1)*sizeof(int));
    LengthPath=(double *)malloc((2*n)*sizeof(double));
    Path=(int *)malloc((2*n)*sizeof(int));

    /* Calcul de l'ordre circulaire X pour D */
    i=1; j=n; odp1(D,X,&i,&j,n);

    /* Initialisation */
    for (i=1; i<=2*n-3; i++)
    {
        MaxCol[i]=0;
        PLACE[i]=0;
        for (j=1;j<=n;j++)
            B[i][j]=0;
    }
    B[1][X[2]]=1; MaxCol[1]=X[2]; Path[1]=1; PlaceNumber=1;
    PLACE[1]=1; LengthPath[1]=D[X[1]][X[2]]; EdgeNumberPath=1; m=1;

    /* Boucle principale de l'algorithme (voir Makarenkov et al.) */
    for(k=2;k<=n-1;k++)
    {
        /* Point 2.1 de l'algorithme */
        DIS=(D[X[1]][X[k]]+D[X[k]][X[k+1]]-D[X[1]][X[k+1]])/2;
        DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;

        if ((DIS<=-EPS1)||(DIS1<=-EPS1)) { printf("\n This is not an additive distance \n");
        free(MaxCol);free(X);free(LengthPath);free(Path);exit(1);return 0; }
        if (DIS<=EPS) {
            DIS=0.0;
        }
        if (DIS1<=EPS) {
            DIS1=0.0;
        }

        S=0.0; i=EdgeNumberPath; if (LengthPath[i]==0.0) i--;
        while (S<=DIS-EPS)
        {
            if (i==0) { S=DIS; break; }
            S=S+LengthPath[i];
            i--;
        }

        /* Point 2.2 de l'algorithme */
        if (fabs(S-DIS)<=EPS)
        {
            M=m+2; DIS=S;
            if (i==0) F=1;
            else if (i==EdgeNumberPath) F=2;
            else { M--; F=3; }
        }
        else {M=m+2; F=0;}

        if (M==m+2)
        {
            if (F==0) { uv=Path[i+1]; EdgeNumberPath=i+2; LengthPath[i+1]=S-DIS; LengthPath[i+2]=DIS1;
            Path[i+1]=m+2; Path[i+2]=m+1;}
            else if (F==1) { uv=Path[1]; EdgeNumberPath=2; LengthPath[1]=0.0; LengthPath[2]=DIS1;
            Path[1]=m+2; Path[2]=m+1;}
            else if (F==2) { uv=Path[EdgeNumberPath]; EdgeNumberPath=EdgeNumberPath+1;LengthPath[EdgeNumberPath]=DIS1;
            Path[EdgeNumberPath-1]=m+2; Path[EdgeNumberPath]=m+1; }

            for (j=1;j<=n;j++)
                B[m+2][j]=B[uv][j];
            MaxCol[m+2]=MaxCol[uv];
        }

        else
        {
            EdgeNumberPath=i+1; LengthPath[i+1]=DIS1; Path[i+1]=m+1;
        }

        /* Point 2.3 de l'algorithme */
        for (j=1;j<=EdgeNumberPath;j++)
            B[Path[j]][X[k+1]]=1;

        /* Point 2.4 de l'algorithme */
        for (j=1;j<=EdgeNumberPath;j++)
            if (MaxCol[Path[j]]<X[k+1]) MaxCol[Path[j]]=X[k+1];

        /* Point 2.5 de l'algorithme */
        for (j=PlaceNumber;j>=1;j--)
            PLACE[j+1]=PLACE[j];
        PLACE[1]=m+1; PlaceNumber++;

        if (M==m+2) {
            i=2;
            while (PLACE[i]!=uv)
                i++;
            for (j=PlaceNumber;j>=i+1;j--)
                PLACE[j+1]=PLACE[j];
            PLACE[i+1]=m+2; PlaceNumber++;
        }

        i=M-1; edge=2;
        do
        {
            if (PLACE[i]==Path[edge])
            {
                edge++; j=i+1;
                while (X[k+1]>MaxCol[PLACE[j]])
                    j++;
                if (j>i+1)
                {
                    l1=PLACE[i];
                    for (l=i+1;l<=j-1;l++)
                        PLACE[l-1]=PLACE[l];
                    PLACE[j-1]=l1;

                }
            }
            i--;
        } while (i!=0);

        m=M;
    }

    /* Libération mémoire */
    free(MaxCol);
    free(X);
    free(LengthPath);
    free(Path);

    return m;
}


/**
 * Compare deux tables de bipartitions (B, B1) et calcule la distance de
 * Robinson-Foulds entre les deux arbres correspondants.
 *
 * La comparaison est effectuée en parcourant les bipartitions triées dans
 * l'ordre de PLACE et PLACE1. Une bipartition commune aux deux arbres
 * incrémente le compteur RF (bipartitions partagées). La distance retournée
 * est la somme des bipartitions non partagées dans chaque arbre.
 *
 * Paramètres :
 *   B, B1         - tables de bipartitions des deux arbres
 *   PLACE, PLACE1 - ordres de tri associés
 *   m, m1         - nombres de bipartitions dans chaque table
 *   n             - nombre de feuilles
 *
 * Retourne la distance de Robinson-Foulds (entier).
 */
int Table_Comparaison (int **B, int ** B1, int *PLACE, int *PLACE1, int m, int m1,int n)
{
    int RF=0,i,p,p1;

    p=1; p1=1;

    while ((p<=m)&&(p1<=m1))
    {
        i=n;
        while ((B[PLACE[p]][i]==B1[PLACE1[p1]][i])&&(i>1))
            i--;
        if (i==1) { RF=RF+1; p++; p1++; }
        else if (B[PLACE[p]][i]>B1[PLACE1[p1]][i]) p1++;
        else p++;

    }

    RF=(m-RF)+(m1-RF);

    return RF;
}


/**
 * Cherche et retourne le premier sommet adjacent à `sommet` dans la matrice
 * d'adjacence dont la distance est inférieure à INFINI.
 *
 * L'arête trouvée est marquée comme visitée en la mettant à INFINI dans les
 * deux sens pour éviter de la revisiter lors des appels récursifs.
 *
 * Paramètres :
 *   Adjacence - matrice d'adjacence (modifiée)
 *   sommet    - indice du sommet courant
 *   n         - nombre d'espèces
 *
 * Retourne l'indice du fils trouvé, ou -1 si aucun n'est accessible.
 */
int findFils(double ** Adjacence,int sommet,int n){

    int i;

    for(i=1;i<=2*n-2;i++){
        if(Adjacence[sommet][i]<INFINI){
            Adjacence[sommet][i] = Adjacence[i][sommet] = INFINI;
            return i;
        }
    }
    return -1;
}


/**
 * Construit récursivement un sous-arbre enraciné en `sommet` à partir de la
 * matrice d'adjacence et du tableau d'arêtes ARETE.
 *
 * Pour chaque sommet, la fonction cherche tous ses fils accessibles via
 * findFils(), puis se rappelle récursivement pour construire leurs sous-arbres.
 * La mémoire allouée doit être libérée via viderArbre().
 *
 * Paramètres :
 *   ARETE     - tableau des arêtes de l'arbre
 *   indice    - compteur de noeuds (partagé entre appels récursifs)
 *   Adjacence - matrice d'adjacence (modifiée lors du parcours)
 *   sommet    - sommet racine du sous-arbre à construire
 *   n         - nombre d'espèces
 *
 * Retourne un pointeur vers la structure TNoeud allouée, ou NULL si sommet == -1.
 */
struct TNoeud * CreerSousArbre(long int *ARETE,int *indice, double ** Adjacence,int sommet,int n){

    int i=0;
    int *tableau;
    int nbElt=0;
    int nouveauFils;
    struct TNoeud * node;

    if(sommet==-1)
        return NULL;

    tableau = safe_malloc<int>(100);

    node = safe_malloc<struct TNoeud>(1);
    node->NoNoeud = sommet;
    node->nbfils = 0;
    node->fils = safe_malloc<struct TNoeud*>(100);
    for(i=0;i<100;i++)
        node->fils[i] = NULL;

    nouveauFils = -1;

    do{
        nouveauFils = findFils(Adjacence,sommet,n);
        if(nouveauFils != -1){
            tableau[nbElt] = nouveauFils;
            nbElt++;
        }
    }while(nouveauFils != -1);

    for(i=0;i<nbElt;i++){
        node->fils[node->nbfils] = CreerSousArbre(ARETE,indice,Adjacence,tableau[i],n);
        node->nbfils = node->nbfils + 1;
    }
    free(tableau);
    return node;
}


/**
 * Trie par ordre croissant les éléments du tableau tab entre les indices
 * debut et fin inclus, par algorithme d'échange simple (tri à bulles).
 *
 * Paramètres :
 *   tab   - tableau d'entiers à trier (modifié en place)
 *   debut - indice de début (inclus)
 *   fin   - indice de fin (inclus)
 */
void sortIntTab(int * tab, int debut, int fin){

    int i,j,tmp;

    for(i=debut;i<=fin;i++){
        for(j=i+1;j<=fin;j++){
            if(tab[i] > tab[j]){
                tmp = tab[i];
                tab[i] = tab[j];
                tab[j] = tmp;
            }
        }
    }
}


/**
 * Parcourt récursivement l'arborescence unNoeud pour remplir la table
 * SousMatriceTree avec les ensembles de sommets de chaque sous-arbre.
 *
 * Pour chaque nœud interne, la fonction agrège les tableaux de sommets de
 * ses fils, les trie, et stocke le résultat dans SousMatriceTree[noNoeud].
 * Pour une feuille, seul son propre numéro est enregistré.
 *
 * Paramètres :
 *   unNoeud        - racine du sous-arbre courant
 *   SousMatriceTree - table à remplir (une entrée par sommet)
 */
void ParcoursArbre(struct TNoeud * unNoeud, struct DescTree * SousMatriceTree){

    int j,k;
    int somme=0;
    int nbSommet=0;
    int nbfils;

    if(unNoeud != NULL){
        nbfils = unNoeud->nbfils;

        if(nbfils != 0){

            for(j=0; j<nbfils;j++){
                ParcoursArbre(unNoeud->fils[j],SousMatriceTree);
            }

            for(j=0; j<nbfils;j++){
                somme = somme + SousMatriceTree[unNoeud->fils[j]->NoNoeud].nbSommet;
            }
            SousMatriceTree[unNoeud->NoNoeud].nbSommet = somme;
            SousMatriceTree[unNoeud->NoNoeud].Tableau = safe_malloc<int>(somme+1);

            for(j=0;j<nbfils;j++){
                for(k=1;k<=SousMatriceTree[unNoeud->fils[j]->NoNoeud].nbSommet;k++) {
                    nbSommet++;
                    SousMatriceTree[unNoeud->NoNoeud].Tableau[nbSommet] = SousMatriceTree[unNoeud->fils[j]->NoNoeud].Tableau[k];
                }
            }
            sortIntTab(SousMatriceTree[unNoeud->NoNoeud].Tableau,1,nbSommet);
            SousMatriceTree[unNoeud->NoNoeud].nbSommet = nbSommet;

        }
        else{
            SousMatriceTree[unNoeud->NoNoeud].Tableau = safe_malloc<int>(2);
            SousMatriceTree[unNoeud->NoNoeud].Tableau[1] = unNoeud->NoNoeud;
            SousMatriceTree[unNoeud->NoNoeud].nbSommet = 1;
        }
    }
}


/**
 * Libère récursivement toute la mémoire allouée pour l'arbre pointé par A.
 *
 * Les enfants sont détruits avant le nœud courant (parcours post-ordre).
 * Le tableau de fils et le nœud lui-même sont libérés à chaque niveau.
 *
 * Paramètre :
 *   A - racine de l'arbre à libérer
 */
void viderArbre(struct TNoeud * A){
    int i;

    if(A->nbfils != 0){
        for(i=0; i<A->nbfils;i++){
            viderArbre(A->fils[i]);
        }
        free(A->fils);
        free(A);
    }
    else{
        free(A->fils);
        free(A);
    }
}


/**
 * Affiche l'arbre enraciné en A sur la sortie standard en représentation
 * indentée (format console).
 *
 * La profondeur prof contrôle le niveau d'indentation. Les sous-arbres
 * sont affichés alternativement avant et après le nœud courant pour
 * obtenir une vue triangulaire équilibrée.
 *
 * Paramètres :
 *   A    - racine de l'arbre à afficher
 *   prof - profondeur courante (0 pour la racine)
 */
void AfficherArbre(struct TNoeud * A,int prof){

    int i, taille;

    if(A != NULL){
        if(A->nbfils== 0){
            for(i=0;i<prof;i++) printf("  ");
            printf("-- %d\n",A->NoNoeud);
        }
        else{
            taille = A->nbfils;
            for(i=0; i<taille/2;i++){
                AfficherArbre(A->fils[i],prof+2);
                printf("\n");
            }
            for(i=0;i<prof;i++) printf("  ");
            printf("-- %d\n",A->NoNoeud);
            for(i=taille/2; i<taille;i++){
                AfficherArbre(A->fils[i],prof+2);
            }
        }
    }
}


/**
 * Convertit un entier non signé val en chaîne de caractères selon la base radix.
 *
 * Si is_neg est vrai, un signe '-' est ajouté en tête. Le résultat est placé
 * dans buf (fourni par l'appelant, doit être assez grand). Les chiffres sont
 * d'abord générés à l'envers puis inversés.
 *
 * Fonction interne utilisée uniquement par itoa_().
 */
static void xtoa (unsigned long val,char *buf,unsigned radix,int is_neg){
    char *p;
    char *firstdig;
    char temp;
    unsigned digval;

    p = buf;

    if (is_neg) {
        /* valeur négative : ajouter le signe et convertir en positif */
        *p++ = '-';
        val = (unsigned long)(-(long)val);
    }

    firstdig = p;

    do {
        digval = (unsigned) (val % radix);
        val /= radix;

        /* conversion en caractère ASCII */
        if (digval > 9)
            *p++ = (char) (digval - 10 + 'a');
        else
            *p++ = (char) (digval + '0');
    } while (val > 0);

    /* Les chiffres sont dans l'ordre inverse : on les retourne */
    *p-- = '\0';

    do {
        temp = *p;
        *p = *firstdig;
        *firstdig = temp;
        --p;
        ++firstdig;
    } while (firstdig < p);
}


/**
 * Convertit l'entier val en chaîne de caractères selon la base radix.
 *
 * Wrapper autour de xtoa() gérant le signe pour la base 10.
 * Le résultat est écrit dans buf (doit être assez grand pour contenir
 * le texte et le '\0' final).
 *
 * Retourne un pointeur vers buf.
 */
char * itoa_(int val,char *buf,int radix){
    if (radix == 10 && val < 0)
        xtoa((unsigned long)val, buf, radix, 1);
    else
        xtoa((unsigned long)(unsigned int)val, buf, radix, 0);
    return buf;
}


/**
 * Filtre les deux matrices de distances pour ne conserver que les espèces
 * présentes dans les deux arbres simultanément.
 *
 * Pour chaque espèce absente de l'un des deux arbres, toutes ses distances
 * sont mises à -1 et son nom est vidé (""). Cette fonction modifie les deux
 * matrices et leurs tableaux de noms en place.
 *
 * Paramètres :
 *   dissSpecies  - matrice de distances de l'arbre espèce (modifiée)
 *   dissGene     - matrice de distances de l'arbre gène (modifiée)
 *   nomsSpecies  - noms d'espèces de l'arbre espèce (modifié)
 *   nomsGene     - noms d'espèces de l'arbre gène (modifié)
 *   nbSpecies    - nombre d'espèces dans l'arbre espèce
 *   nbGene       - nombre d'espèces dans l'arbre gène
 *   fichier      - chemin du fichier d'erreurs (non utilisé ici)
 */
void filtrerMatrice(double **dissSpecies, double **dissGene, char **nomsSpecies, char **nomsGene,int nbSpecies, int nbGene, char * fichier){

        int i,j,temoin;

        /* Supprimer de l'arbre espèce les espèces absentes de l'arbre gène */
        for(i = 1 ; i <= nbSpecies; i++){
            temoin = 0;
            for(j=1;j<=nbGene;j++){
                if(strcmp(nomsSpecies[i],nomsGene[j])==0)
                    temoin=1;
            }
            if(temoin == 0){
                for(j=1;j<=nbSpecies;j++){
                    dissSpecies[i][j] = dissSpecies[j][i] = -1;
                }
                strcpy(nomsSpecies[i],"");
            }
        }

        /* Supprimer de l'arbre gène les espèces absentes de l'arbre espèce */
        for(i = 1; i <= nbGene; i++){
            temoin = 0;
            for(j=1;j<=nbSpecies;j++){
                if(strcmp(nomsSpecies[j],nomsGene[i])==0)
                    temoin=1;
            }
            if(temoin == 0){
                for(j=1; j <= nbGene; j++){
                    dissGene[i][j] = dissGene[j][i] = -1;
                }
                strcpy(nomsGene[i],"");
            }
        }
}


/**
 * Écrit la matrice de distances mat dans le fichier outfile au format
 * PHYLIP (première ligne : taille, puis une ligne par espèce).
 *
 * Seules les espèces dont le nom n'est pas vide ("") sont écrites.
 * Si le nombre d'espèces communes est inférieur à 3, la fonction retourne
 * -1 sans écrire (arbre non exploitable).
 *
 * Paramètres :
 *   mat     - matrice de distances à écrire
 *   outfile - chemin du fichier de sortie (écrasé si existant)
 *   taille  - nombre total d'espèces dans mat
 *   noms    - tableau des noms d'espèces
 *
 * Retourne le nombre d'espèces écrites, ou -1 si moins de 3 espèces communes.
 */
int ecrireMatrice(double **mat,const char *outfile,int taille,char **noms){
    int i,j,finalTaille;
    FILE *out;

    finalTaille=0;
    for(i=1;i<=taille;i++)
        if(strcmp(noms[i],"")!=0)
            finalTaille = finalTaille+1;
    if(finalTaille < 3){
        return -1;
    }
    out = safe_fopen(outfile, "w+");
    fprintf(out, "%d", finalTaille);
    for (i = 1; i <= taille; i++) {
        if (strcmp(noms[i], "") != 0) {
            fprintf(out, "\n%s", noms[i]);
            for (j = 1; j <= taille; j++)
                if (mat[i][j] != -1) {
                    fprintf(out, " %lf", mat[i][j]);
                }
        }
    }
    fclose(out);
    return finalTaille;
}


/**
 * Ajoute la matrice de distances mat au fichier outfile en mode ajout (append).
 *
 * Une ligne vide est d'abord écrite pour séparer les blocs. Seules les espèces
 * dont le nom n'est pas vide sont ajoutées.
 *
 * Paramètres :
 *   mat     - matrice de distances de l'arbre gène
 *   outfile - chemin du fichier de sortie (ouvert en mode append)
 *   taille  - nombre total d'espèces dans mat
 *   noms    - tableau des noms d'espèces
 */
void ajouterMatriceGene(double **mat, const char *outfile, int taille, char **noms) {
    int i, j;
    FILE *out = safe_fopen(outfile, "a");
    fprintf(out, "\n");
    for (i = 1; i <= taille; i++) {
        if (strcmp(noms[i], "") != 0) {
            fprintf(out, "\n%s", noms[i]);
            for (j = 1; j <= taille; j++)
                if (mat[i][j] != -1)
                    fprintf(out, " %lf", mat[i][j]);
        }
    }
    fclose(out);
}


/**
 * Trie la matrice DISS et son tableau de noms NomsDISS pour que l'ordre
 * des espèces corresponde à celui de NomsADD.
 *
 * Pour chaque espèce dans NomsADD, la fonction cherche sa position dans
 * NomsDISS et réorganise les lignes/colonnes de DISS en conséquence.
 * Utilisée pour aligner les matrices de distances espèce et gène.
 *
 * Paramètres :
 *   DISS      - matrice de distances à réordonner (modifiée)
 *   NomsDISS  - noms associés à DISS (modifiés)
 *   NomsADD   - noms de référence (ordre cible)
 *   n         - nombre d'espèces
 */
void TrierMatrices(double **DISS,char **NomsDISS,char **NomsADD,int n)
{
    int trouve, ligne,colonne,i,j;
    double ** DISS_;
    char   ** NomsDISS_;
    char noms[50];
    int * table;

    table = (int *) malloc((n+1)*sizeof(int));
    DISS_ = (double **) malloc((n+1)*sizeof(double*));
    NomsDISS_ = (char **) malloc((n+1)*sizeof(char*));

    for (i=0;i<=n;i++)
    {
        DISS_[i] = (double*)malloc((n+1)*sizeof(double));
        NomsDISS_[i] = (char*) malloc(SPECIES_NAME_LENGTH);
        if (NomsDISS_[i] == NULL || DISS_[i] == NULL) {
            fprintf(stderr, "memory allocation failed in TrierMatrices\n");
            exit(1);
        }
    }

    for(ligne = 1;ligne<=n;ligne++)
    {
        strcpy(noms,NomsADD[ligne]);
        trouve = 0;
        for(colonne = 1;colonne<=n;colonne++)
        {
            if(strcmp(noms,NomsDISS[colonne])==0)
            {
                trouve = 1;
                table[ligne] = colonne;
                strcpy(NomsDISS_[ligne],noms);
            }
        }
        if(trouve==0)
        {
            printf("\n%s %s",noms,"is not in the gene matrix.This program must stop");
            exit (-1);
        }
    }

    for(i=1;i<=n;i++)
        for(j=1;j<=i;j++) DISS_[i][j] = DISS_[j][i] = DISS[table[i]][table[j]];

    for(i=1;i<=n;i++)
    {
        strcpy(NomsDISS[i],NomsDISS_[i]);
        for(j=1;j<=n;j++) DISS[i][j] = DISS_[i][j];
    }
    for (i=0;i<=n;i++)
    {
        free(DISS_[i]);
        free(NomsDISS_[i]);
    }
    free(DISS_);
    free(NomsDISS_);
    free(table);
}


/**
 * Analyse une chaîne au format Newick et extrait les arêtes (ARETE),
 * leurs longueurs (LONGUEUR) et les noms des feuilles (lesNoms).
 *
 * Le format attendu est : arborescence entre parenthèses, noms de feuilles
 * suivis de ':' et d'une longueur de branche, terminée par ';'.
 *
 * Algorithme : traitement itératif des paires de parenthèses les plus
 * profondes, en construisant progressivement les tableaux ARETE/LONGUEUR.
 * Les noeuds internes sont représentés par des numéros > n (n+1, n+2, ...).
 *
 * Paramètres :
 *   newick  - chaîne Newick de l'arbre
 *   ARETE   - tableau des arêtes résultant (modifié, indexé à partir de 0)
 *   LONGUEUR - longueurs des arêtes (modifié)
 *   lesNoms - noms des feuilles (modifié, indexé à partir de 1)
 *   kt      - nombre de branches terminales (modifié)
 *
 * Retourne le nombre de feuilles n, ou -1 en cas d'erreur de format.
 */
int lectureNewick(string newick, long int * ARETE, double * LONGUEUR, char ** lesNoms, int *kt)
{
    int n;
    int cpt_x;
    int i,j,j1,k, a, a1, a2,a3 = 0, VertexNumber,numero = 0;
    char symbol, *string, *string1, *string2;
    int taxaPos;
    int aretePos;
    char symbolOld =' ';
    int zz, xx,jj,ii;
    double longueur;
    char * tempString;
    int cpt=0;
    int na;
    char *string4 = (char*) malloc(1000);
    int temoin=0;

    a=0;

    /* Vérification du format Newick et comptage des espèces */
    i=0;
    n = 0;

    do
    {
        symbol = newick[cpt++];
        if (symbol==':') i++;
        if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
        if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
        if(symbol==':' && temoin==1) temoin=0;
        symbolOld = symbol;
    }  while(symbol != '\0');

    cpt=0;
    na=i;
    if(i<=2*n-3)(*kt)=i;
    else (*kt)=2*n-3;

    if (i == 0) { printf("Incorrect Newick file format. Edge lengths must be indicated after a ':' characters."); return -1;}

    if ((i>2*n-3) || (i<n)){
        /* Format non standard : continuer sans erreur fatale */
    }

    /* Compter les parenthèses ouvrantes */
    j=0;
    do{
        symbol=newick[cpt++];
        if (symbol=='(') j++;
    }while(symbol != '\0');

    cpt=0;
    j1=0;

    /* Compter les parenthèses fermantes */
    do{
        symbol=newick[cpt++];
        if (symbol==')') j1++;
    }while(symbol != '\0');

    cpt=0;

    /* Vérification de la cohérence des parenthèses */
    if (j1!=j) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); return -1;}

    k=0;

    do{
        symbol=newick[cpt++];
        if (symbol==',') k++;
    }while(symbol != '\0');

    cpt=0;

    a=0;

    do{
        symbol=newick[cpt++];
        if (symbol==';') a++;
    }while(symbol != '\0');
    cpt=0;

    if (a==0) { printf("Incorrect Newick file format. Newick string must be followed by a ';' character."); return -1;}
    else if (a>1) { printf("Incorrect Newick file format. Newick string must contain (in the end) only one ';' character."); return -1;}

    a=0;
    do{
        symbol=newick[cpt++];a++;
    }while(symbol == ' ');

    cpt=0;

    if (symbol!='(') { printf("Incorrect Newick file format. Newick string must begin with a '(' character."); return -1;}

    a=0;

    do{ symbol=newick[cpt++];
        if (symbol=='%') a++;
    }while(symbol != '\0');

    cpt=0;

    if (a>0) { printf("Incorrect Newick file format. Newick string cannot contain '%%' character."); return -1;}

    /* Vérification de la longueur des noms (max 50 caractères) */
    do
    {
        symbol=newick[cpt++];
        if ((symbol=='(')||(symbol==','))
        {
            symbol=newick[cpt++]; a=0;
            if ((symbol!='(')&&(symbol!=',')&&(symbol!=';')&&(symbol!=':'))
            {
                cpt--;
                do{
                    symbol=newick[cpt++];a++;
                }while ((symbol!=':')&&(symbol!='\0'));
            }
            else cpt--;
            if (a>50) { printf("Incorrect Newick file format. Names of objects must not exceed 50 characters.");  return -1;}
        }
    }while(symbol != '\0');
    cpt=0;

    string = (char*)malloc((100000*n) * sizeof(char));
    string2 = (char*)malloc((100000*n) * sizeof(char));
    string1 = (char*)malloc((200000) * sizeof(char));

    if ((string == NULL)||(string1 == NULL)||(string2 == NULL))
    { printf("Input data are too large or not a correct Newick file chosen"); return -1;}

    a=0;

    /* Supprimer les espaces, tabulations et retours à la ligne */
    do{
        symbol=newick[cpt++];
        if ((symbol!=' ')&&(symbol!='\n')&&(symbol!='\t')) { string[a++]=symbol; }
    }while(symbol !='\0');

    int boot_value;
    k=0; VertexNumber=n;
    taxaPos =1;
    aretePos = 1;
    boot_value=0;

    /* Traitement itératif de la chaîne Newick : décomposition des paires () */
    while (string[0] == '(')
    {
        a1 = 0;
        a2 = 0;
        while( string[a2] != ')')
        {
            if(string[a2] == '(') a1 = a2;
            a2++;
        }

        /* a1 : début du nœud à traiter (parenthèse ouvrante)
         * a2 : fin du nœud à traiter (parenthèse fermante)
         * a3 : délimite le nœud et sa longueur de branche */

        zz = a1+1;
        VertexNumber++;
        boot_value=0;
        for ( ii = a1+1; ii <= a2; ii++)
        {
            if (string[ii] == ':')
            {
                xx = 0;
                a3 = ii+1;

                if( string[zz] == '%')
                {
                    /* Nœud interne déjà traité : récupérer son numéro */
                    for ( jj = zz+1; jj < ii; jj++)
                    {
                        if(string[jj] == '|')
                            break;
                        string1[xx++] = string[jj];
                    }
                    string1[xx++] = '\0';
                    numero = atoi(string1);

                    if(string[jj] == '|' ){
                        boot_value=1;
                        cpt_x=0;
                        jj++;
                        while(string[jj] != ':')
                            string4[cpt_x++] = string[jj++];
                        string4[cpt_x] = '\0';
                    }

                }
                else
                {
                    /* Feuille : récupérer le nom du taxon */
                    for(jj = zz; jj < ii; jj++)
                    {
                        lesNoms[taxaPos-1][xx++] = string[jj];
                    }
                    numero = taxaPos;
                    lesNoms[taxaPos-1][xx] = '\0';
                    taxaPos++;
                }

            }
            else if(string[ii] == ','  || string[ii] == ')')
            {
                /* Lire la longueur de branche et enregistrer l'arête */
                xx = 0;
                zz = ii +1;
                for ( jj = a3; jj < ii; jj++)
                {
                    string1[xx++] = string[jj];
                }
                string1[xx++] = '\0';
                longueur = atof(string1);
                ARETE[aretePos++] = VertexNumber;
                ARETE[aretePos++] = numero;
                LONGUEUR[(aretePos-1)/ 2] = longueur;

                boot_value=0;
            }

        }

        /* Transcrire la nouvelle chaîne avec le nœud interne encodé en '%' */
        xx = 0;
        for ( jj = 0; jj < (int)a1; jj++)
        {string2[xx++] = string[jj];}

        itoa_(VertexNumber,string1,10);
        string2[xx++] = '%';
        for( jj = 0; jj < (int) strlen(string1); jj++)
        {string2[xx++] = string1[jj];}

        int temoin=0;

        for( jj = a2+1; jj <= a; jj++)
        {
            if(string[jj] != ':' && temoin==0){
                string2[xx++] = '|';
            }
            temoin = 1;
            string2[xx++] = string[jj];
        }

        /* Échanger string et string2 */
        tempString = string;
        string = string2;
        string2 = tempString;
        tempString = 0;
        a = xx;

    } /* fin du while : toute la chaîne a été traitée */


    int root_existance = -1;

    for(jj=0;jj<n;jj++){
        if (strcmp(lesNoms[jj],"Root") == 0)
            root_existance = jj;
    }

    ARETE[aretePos++] = 0;
    ARETE[aretePos++] = 0;

    /* Décaler ARETE d'une position (alignement base 0 → base 1) */
    i=0;
    int cpt_branches=0;
    do{
        i++;cpt_branches++;
        ARETE[i-1] = ARETE[i];
        i++;
        ARETE[i-1] = ARETE[i];
    }while(ARETE[i] != 0);

    for(i=1;i<=na;i++){
        LONGUEUR[i-1] = LONGUEUR[i];
    }

    /* Gestion de la racine explicite ("Root") */
    if( root_existance > 0){
        long noeud_interne = -1.0;
        for(i=1;i<=na;i++){
            if(ARETE[2*i-1] == root_existance)
                noeud_interne = ARETE[2*i-2];
            if(ARETE[2*i-2] == root_existance)
                noeud_interne = ARETE[2*i-1];
        }

        for(i=1;i<=na;i++){
            if((ARETE[2*i-1] != root_existance) && (noeud_interne == ARETE[2*i-2])){
                LONGUEUR[i-1] = 50;
            }
            if((ARETE[2*i-2] != root_existance) && (noeud_interne == ARETE[2*i-1])){
                LONGUEUR[i-1] = 50;
            }
        }
    }

    /* Détection et suppression du nœud de degré 2 (pseudo-racine) */
    int * tableau = (int*)malloc((2*n)*sizeof(int));
    int deg2=-1,deg1=-1;
    for(i=1;i<=2*n;i++)
        tableau[i-1] = 0;
    for(i=1;i<=na;i++){
        tableau[ARETE[2*i-1]]++;
        tableau[ARETE[2*i-2]]++;
    }

    i=n+1;
    while(tableau[i] > 0){
        if(tableau[i] == 2) deg2 = i;
        if(tableau[i] == 1) deg1 = i;
        i++;
    }
    int pos_racine=-1;

    for(i=1;i<=na;i++){
        if(ARETE[2*i-1] == deg1){ ARETE[2*i-1] = deg2; pos_racine=i;  break;}
        if(ARETE[2*i-2] == deg1){ ARETE[2*i-2] = deg2; pos_racine=i;  break;}
    }
    if(pos_racine != -1){
        LONGUEUR[pos_racine-1] = 100;
    }
    else if(deg2 != -1){
        int pos1=-1,pos2=-1;
        for(i=1;i<=na;i++){
            if((ARETE[2*i-1] == deg2 || ARETE[2*i-2] == deg2) && pos1==-1){ pos1=i;}
            else if((ARETE[2*i-1] == deg2 || ARETE[2*i-2] == deg2) && pos1 != -1){ pos2=i;}
        }

        if(ARETE[2*pos1-1] == deg2)
            ARETE[2*pos1-1] = (ARETE[2*pos2-1] == deg2)?ARETE[2*pos2-2]:ARETE[2*pos2-1];
        else
            ARETE[2*pos1-2] = (ARETE[2*pos2-1] == deg2)?ARETE[2*pos2-2]:ARETE[2*pos2-1];

        for(i=pos2;i<=na;i++){
            LONGUEUR[i-1] = LONGUEUR[i];
            ARETE[2*i-1] = ARETE[2*(i+1)-1];
            ARETE[2*i-2] = ARETE[2*(i+1)-2];
        }
        na--;
    }

    (*kt) = 2*n-3 - (*kt);

    /* Décaler les noms de feuilles : lesNoms[0] → lesNoms[1] */
    for(i = n; i >= 1; i--){
        strcpy(lesNoms[i],lesNoms[i-1]);
    }
    return n;
}


/* =========================================================================
 * Fonctions utilitaires de débogage — par Arthur Debeaupte
 * ========================================================================= */

/**
 * Écrit le contenu d'une structure InputTree dans un fichier texte situé
 * dans le répertoire "../treeFiles/".
 *
 * Tous les champs de la structure sont écrits (size, ADD, Input, SpeciesName,
 * Root, Adjacence, ARETE, LONGUEUR, W, kt, degre). Les tableaux NULL sont
 * signalés explicitement.
 *
 * Paramètres :
 *   filename - nom du fichier de sortie (doit se terminer par ".txt")
 *   tree     - structure InputTree à écrire
 */
void writeInputTreeToFile(const std::string& filename, InputTree& tree) {
    std::string finalFileName = "../treeFiles/" + filename;
    std::ofstream file(finalFileName);

    if (!file) {
        std::cerr << "Failed to create the file: " << filename << std::endl;
        return;
    }

    file << "size : " << tree.size << std::endl;

    file << "\nADD : " << std::endl;
    if (tree.ADD != NULL) {
        for (int i = 0; i < tree.size; i++) {
            for (int j = 0; j < tree.size; j++) {
                file << tree.ADD[i][j] << " ";
            }
            file << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "\nInput : " << std::endl;
    if (tree.Input != NULL) {
        for (int i = 0; i < tree.size; i++) {
            for (int j = 0; j < tree.size; j++) {
                file << tree.Input[i][j] << " ";
            }
            file << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "\nSpeciesName : " << std::endl;
    if (tree.SpeciesName != NULL) {
        for (int i = 0; i < tree.size; i++) {
            file << tree.SpeciesName[i] << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "\nRoot : " << tree.Root << std::endl;

    file << "\nAdjacence : ";
    if (tree.Adjacence != NULL ) {
        for (int i = 0; i < tree.size; i++) {
            for (int j = 0; j < tree.size; j++) {
                file << std::fixed << std::setprecision(2) << tree.Adjacence[i][j] << " ";
            }
            file << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "\nARETE : ";
    if (tree.ARETE != NULL) {
        for (int i = 0; i < tree.size; i++) {
            file << tree.ARETE[i] << " ";
        }
        file << std::endl;
    } else {
        file << " NULL" << std::endl;
    }

    file << "\nLONGUEUR : ";
    if (tree.LONGUEUR != NULL) {
        for (int i = 0; i < tree.size; i++) {
            file << tree.LONGUEUR[i] << " ";
        }
        file << std::endl;
    } else {
        file << "NULL" << std::endl;
    }

    file << "\nW :" << std::endl;
    if (tree.W != NULL) {
        for (int i = 0; i < tree.size; i++) {
            for (int j = 0; j < tree.size; j++) {
                file << tree.W[i][j] << " ";
            }
            file << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "\nkt : " << tree.kt << std::endl;

    file << "\ndegre : ";
    if (tree.degre != NULL) {
        for (int i = 0; i < tree.size; i++) {
            file << tree.degre[i] << " ";
        }
        file << std::endl;
    } else {
        file << "NULL" << std::endl;
    }

    file.close();
}


/**
 * Écrit le contenu d'une structure CRITERIA dans un fichier texte situé
 * dans le répertoire "../treeFiles/".
 *
 * Les champs scalaires (LS, rLS, BD, rBD, rRF, diff_bd, RF, QD, m,
 * nbHgtFound) et les matrices B, BI, PLACE sont écrits. Les tableaux NULL
 * sont signalés explicitement.
 *
 * Paramètres :
 *   filename - nom du fichier de sortie (doit se terminer par ".txt")
 *   criteria - structure CRITERIA à écrire
 */
void writeCriteriaToFile(const std::string& filename, CRITERIA& criteria) {
    std::string finalFileName = "../treeFiles/" + filename;
    std::ofstream file(finalFileName);

    if (!file) {
        std::cerr << "Failed to create the file: " << filename << std::endl;
        return;
    }

    file << "LS : " << criteria.LS << std::endl;
    file << "rLS : " << criteria.rLS << std::endl;
    file << "BD : " << criteria.BD << std::endl;
    file << "rBD : " << criteria.rBD << std::endl;
    file << "rRF : " << criteria.rRF << std::endl;
    file << "diff_bd : " << criteria.diff_bd << std::endl;
    file << "RF : " << criteria.RF << std::endl;
    file << "QD : " << criteria.QD << std::endl;
    file << "m : " << criteria.m << std::endl;
    file << "nbHgtFound : " << criteria.nbHgtFound << std::endl;

    file << "B : " << std::endl;
    if (criteria.B != NULL) {
        for (int i = 0; i < criteria.m; i++) {
            for (int j = 0; j < criteria.m; j++) {
                file << criteria.B[i][j] << " ";
            }
            file << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "BI : " << std::endl;
    if (criteria.BI != NULL) {
        for (int i = 0; i < criteria.m; i++) {
            for (int j = 0; j < criteria.m; j++) {
                file << criteria.BI[i][j] << " ";
            }
            file << std::endl;
        }
    } else {
        file << "NULL" << std::endl;
    }

    file << "PLACE : ";
    if (criteria.PLACE != NULL) {
        for (int i = 0; i < criteria.m; i++) {
            file << criteria.PLACE[i] << " ";
        }
        file << std::endl;
    } else {
        file << "NULL" << std::endl;
    }

    file.close();
}
