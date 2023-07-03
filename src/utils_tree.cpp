//
//  utils_tree.cpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

#include "utils_tree.hpp"


//==============================================
//=
//==============================================

double BipartitionDistance (int **B, int ** B1,int n)
{
    int i,j,k,k1,cpt1,cpt2,cpt3,cpt4,t1,t2;
    double c1,c2,min1 = 0.0,min2 = 0.0,cptB=0.0,cptB1=0.0;
    int *flag,*flag1,**Bi,**Bi1;
    Bi=(int **) malloc((2*n-2)*sizeof(int*));
    Bi1=(int **) malloc((2*n-2)*sizeof(int*));
    for (i=0;i<2*n-2;i++)
    {
        Bi[i]=(int *) malloc((2*n)*sizeof(int));
        Bi1[i]=(int *) malloc((2*n)*sizeof(int));
        if ((Bi1[i]==NULL)||(Bi[i]==NULL))
        {
            printf(" Data matrix is too large\n ");
            exit(2);
        }
    }
    flag = (int*) malloc(2*n*(sizeof(int)));
    flag1 = (int*) malloc(2*n*(sizeof(int)));
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
    
    for(i=1;i<=n-3;i++){
        t1=t2=0;
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


//===========================================================
//== Floyd
//===========================================================

void Floyd(double ** Adjacence , double ** DIST,int n,int kt)
{
    int i,j,k;

    for(i=1;i<=2*n-2-kt;i++)
        for(j=1;j<=2*n-2-kt;j++)
        {
            if(i==j)
                DIST[i][j] = 0;
            else
                DIST[i][j] = Adjacence[i][j];
        }

        for(i=1;i<=2*n-2-kt;i++)
            for(j=1;j<=2*n-2-kt;j++)
                for(k=1;k<=2*n-2-kt;k++)
                {
                    if((DIST[j][i] + DIST[i][k]) < DIST[j][k])
                        DIST[j][k] = DIST[j][i] + DIST[i][k];
                }
}

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

        for(i=1;i<=2*n-2-kt;i++)
            for(j=1;j<=2*n-2-kt;j++)
                for(k1=1;k1<=2*n-2-kt;k1++)
                {
                    if((DIST[j][i] + DIST[i][k1]) < DIST[j][k1])
                        DIST[j][k1] = DIST2[j][k1] = DIST[j][i] + DIST[i][k1];
                }
}


//===========================================================
//==
//===========================================================


void loadAdjacenceMatrix( double **Adjacence, long int *ARETE, double *LONGUEUR,int size,int kt){
    
    int i,j;
    
    for(i=1;i<=2*size-2;i++) /*/(n+1)*/
        for(j=1;j<=2*size-2;j++){
            Adjacence[i][j] = Adjacence[j][i] = INFINI;
}
    for(i=1;i<=2*size-3-kt;i++){
        Adjacence[ARETE[2*i-2]][ARETE[2*i-1]] = LONGUEUR[i-1];//(LONGUEUR[i-1]>5*epsilon)?LONGUEUR[i-1]:5*epsilon;
        Adjacence[ARETE[2*i-1]][ARETE[2*i-2]] = LONGUEUR[i-1]; //(LONGUEUR[i-1]>5*epsilon)?LONGUEUR[i-1]:5*epsilon;
    }
}

//=====================================================
//==
//=====================================================

void odp1(double **D, int *X, int *i1, int *j1, int n)
{
    double S1,S;
    int i,j,k = 0,a,*Y1;

    Y1=(int *) malloc((n+1)*sizeof(int));

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


//=====================================================
//==
//=====================================================

int Tree_edges (double **DI, long int *ARETE, double *LONGUEUR, int n,int binaire)
{

    struct EDGE { unsigned int U; unsigned int V; double LN;};
    struct EDGE *Path,*Tree;
    int i,j,k,p,P,*X;
    double S,DIS,DIS1,*L,**D;
    int pasfini = 1;
    int kt=0;
    long SomToDel,OtherSom;

    X=(int *)malloc((n+1)*sizeof(int));
    L=(double *)malloc((n+1)*sizeof(double));
    Tree=(struct EDGE *)malloc((2*n-2)*sizeof(struct EDGE));
    Path=(struct EDGE *)malloc((n+2)*sizeof(struct EDGE));


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

    /* Verification de l'equivalence des topologies */
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
    /*/== rajoute le 22 avril 2005*/
    while((pasfini==1)&&(binaire==0)){
        pasfini = 0;
        for (i=1;i<=2*n-3-kt;i++){
            if((LONGUEUR[i-1] == 2*epsilon_ba)&&(ARETE[2*i-2]>n)&&(ARETE[2*i-1]>n)){    //= branche interne de taille=2*epsilon
                if(ARETE[2*i-2] > ARETE[2*i-1]){
                    SomToDel = ARETE[2*i-2];
                    OtherSom = ARETE[2*i-1];
                }
                else{
                    SomToDel = ARETE[2*i-1];
                    OtherSom = ARETE[2*i-2];
                }
                pasfini=1;

                //== trouver les branches connexes*/
                for (j=1;j<=2*n-3-kt;j++){

                    if(j!=i){
                        if(ARETE[2*j-2] == SomToDel){/*&&(ARETE[2*j-1] != ARETE[2*i-1])){*/
                            ARETE[2*j-2] = OtherSom;
                        }
                        else if(ARETE[2*j-1] == SomToDel){/*/&&(ARETE[2*j-2] != ARETE[2*i-1])){*/
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


//==============================================================
//
//==============================================================


int Bipartition_Table (double **D, int **B, int *PLACE, int n)
{
    int i,j,k,l,l1,*MaxCol,*X,EdgeNumberPath,m,uv = 0,PlaceNumber,edge,*Path,M,F;
    double S,DIS,DIS1,*LengthPath;
    double EPS=1.e-5;
    double EPS1=1.e-2;

    /* Memory allocation */

    MaxCol=(int *)malloc((2*n)*sizeof(int));
    X=(int *)malloc((2*n+1)*sizeof(int));
    LengthPath=(double *)malloc((2*n)*sizeof(double));
    Path=(int *)malloc((2*n)*sizeof(int));

    /* Computation of a circular order X for D */

    i=1; j=n; odp1(D,X,&i,&j,n);

    /* Initialization */
    for (i=1; i<=2*n-3; i++)
    {
        MaxCol[i]=0;
        PLACE[i]=0;
        for (j=1;j<=n;j++)
            B[i][j]=0;
    }
    B[1][X[2]]=1; MaxCol[1]=X[2]; Path[1]=1; PlaceNumber=1;
    PLACE[1]=1; LengthPath[1]=D[X[1]][X[2]]; EdgeNumberPath=1; m=1;

    /* The main loop */

    for(k=2;k<=n-1;k++)
    {
        /* Point 2.1 of the algorithm (see the referenced article by Makarenkov and Leclerc) */

        DIS=(D[X[1]][X[k]]+D[X[k]][X[k+1]]-D[X[1]][X[k+1]])/2;
        DIS1=(D[X[1]][X[k+1]]+D[X[k]][X[k+1]]-D[X[1]][X[k]])/2;
        
        if ((DIS<=-EPS1)||(DIS1<=-EPS1)) { printf("\n This is not an additive distance \n");
        free(MaxCol);free(X);free(LengthPath);free(Path);exit(1);return 0; }
        if (DIS<=EPS) DIS=0.0; if (DIS1<=EPS) DIS1=0.0;

        S=0.0; i=EdgeNumberPath; if (LengthPath[i]==0.0) i--;
        while (S<=DIS-EPS)
        {
            if (i==0) { S=DIS; break; }  /* checking the limit */
            S=S+LengthPath[i];
            i--;
        }

        /* Point 2.2 of the algorithm */

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

        /* Point 2.3 of the algorithm */

        for (j=1;j<=EdgeNumberPath;j++)
            B[Path[j]][X[k+1]]=1;

        /* Point 2.4 of the algorithm */

        for (j=1;j<=EdgeNumberPath;j++)
            if (MaxCol[Path[j]]<X[k+1]) MaxCol[Path[j]]=X[k+1];

        /* Point 2.5 of the algorithm */

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


    /* memeory liberation */
    free(MaxCol);
    free(X);
    free(LengthPath);
    free(Path);

    return m;

}

//==============================================================
//
//==============================================================


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


//=======================================
//
//=======================================

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

//=======================================
//
//=======================================
 
struct TNoeud * CreerSousArbre(long int *ARETE,int *indice, double ** Adjacence,int sommet,int n){

    int i=0;
    int *tableau;
    int nbElt=0;
    int nouveauFils;
    struct TNoeud * node;

    if(sommet==-1)
        return NULL;

    tableau = (int*)malloc(100*sizeof(int));

    node = (struct TNoeud*)malloc(sizeof(struct TNoeud));
    node->NoNoeud = sommet;
    node->nbfils = 0;
    node->fils = (struct TNoeud**)malloc(100*sizeof(struct TNoeud*));
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


//=============================================
//= tri bulle d'un tableau d'entiers
//=============================================

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


//============================================================
//
//============================================================


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
            SousMatriceTree[unNoeud->NoNoeud].nbSommet = somme; /*SousMatriceTree[unNoeud->droit->NoNoeud].nbSommet + SousMatriceTree[unNoeud->gauche->NoNoeud].nbSommet;*/
            SousMatriceTree[unNoeud->NoNoeud].Tableau = (int*)malloc((somme+1)*sizeof(int));
            
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

            SousMatriceTree[unNoeud->NoNoeud].Tableau = (int*)malloc((2)*sizeof(int));

            SousMatriceTree[unNoeud->NoNoeud].Tableau[1] = unNoeud->NoNoeud;
            SousMatriceTree[unNoeud->NoNoeud].nbSommet = 1;
        }
    }
}


//=======================================
//
//=======================================


void viderArbre(struct TNoeud * A){
    int i;
    
    if(A->nbfils != 0){
        for(i=0; i<A->nbfils;i++){
            viderArbre(A->fils[i]);
        }
    //    printf("%d ",A->NoNoeud);
        free(A->fils);
        free(A);
    }
    else{
        free(A->fils);
        free(A);
    }
}


//=======================================
//
//=======================================


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


//=======================================
//
//=======================================


static void xtoa (unsigned long val,char *buf,unsigned radix,int is_neg){
    char *p;                /* pointer to traverse string */
    char *firstdig;         /* pointer to first digit */
    char temp;              /* temp char */
    unsigned digval;        /* value of digit */

    p = buf;

    if (is_neg) {
        /* negative, so output '-' and negate */
        *p++ = '-';
        val = (unsigned long)(-(long)val);
    }

    firstdig = p;           /* save pointer to first digit */

    do {
        digval = (unsigned) (val % radix);
        val /= radix;       /* get next digit */

        /* convert to ascii and store */
        if (digval > 9)
            *p++ = (char) (digval - 10 + 'a');  /* a letter */
        else
            *p++ = (char) (digval + '0');       /* a digit */
    } while (val > 0);

    /* We now have the digit of the number in the buffer, but in reverse
    order.  Thus we reverse them now. */

    *p-- = '\0';            /* terminate string; p points to last digit */

    do {
        temp = *p;
        *p = *firstdig;
        *firstdig = temp;   /* swap *p and *firstdig */
        --p;
        ++firstdig;         /* advance to next two digits */
    } while (firstdig < p); /* repeat until halfway */
}


/* Actual functions just call conversion helper with neg flag set correctly,
and return pointer to buffer. */

char * itoa_(int val,char *buf,int radix){
    if (radix == 10 && val < 0)
        xtoa((unsigned long)val, buf, radix, 1);
    else
        xtoa((unsigned long)(unsigned int)val, buf, radix, 0);
    return buf;
}



//=======================================
//
//=======================================


void filtrerMatrice(double **dissSpecies, double **dissGene, char **nomsSpecies, char **nomsGene,int nbSpecies, int nbGene, char * fichier){

        int i,j,temoin;
        
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


//== write the matrix in the output file

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
    if((out=fopen(outfile,"w+"))==NULL){
        printf("cannot create output file !!");
        printf("==1==>%s<===",outfile);
        exit(-1);
    }
    else{
        fprintf(out,"%d",finalTaille);
        for(i=1;i<=taille;i++){
            if(strcmp(noms[i],"") != 0){//if(strlen(noms[i]) > 1){
                fprintf(out,"\n%s",noms[i]);
                for(j=1;j<=taille;j++)
                    if(mat[i][j] != -1){
                        fprintf(out," %lf",mat[i][j]);
                    }
            }
        }
        fclose(out);
    }
    return finalTaille;
}


//== add the gene matrix to the input file

void ajouterMatriceGene(double **mat,const char *outfile,int taille,char **noms) {
    int i,j;
    FILE *out;
    if((out=fopen(outfile,"a"))==NULL){
        printf("cannot create output file !!");
        printf("==2==>%s<===",outfile);
        exit(-1);
    }
    else{
        fprintf(out,"\n");
        for(i=1;i<=taille;i++){
            if(strcmp(noms[i],"") != 0){  //if(strlen(noms[i]) > 1){
                fprintf(out,"\n%s",noms[i]);
                for(j=1;j<=taille;j++)
                    if(mat[i][j] != -1)
                        fprintf(out," %lf",mat[i][j]);
            }
        }
        fclose(out);
    }
}


/*
* This procedure sorts the 2nd matrix read according to
* of the name of the species of the first matrix
* DISS: matrix of the gene
* NomsDISS: array of species names of gene matrix
* NomsADD: table containing the species names of the species matrix
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
        DISS_[i]=(double*)malloc((n+1)*sizeof(double));
        NomsDISS_[i] = (char*) malloc((n+1)*sizeof(15));
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


//===========================================================================
//const char * newick
//===========================================================================
//===========================================================================


int lectureNewick(string newick, long int * ARETE, double * LONGUEUR, char ** lesNoms, int *kt)
{
    int n;
    int cpt_x;
    // Ce sous programme permet de lire un arbre au format newick et de le transcrire dans
    // des vecteurs arete-longueur commencant à 1
    // ATTENTION: les noms commencent à 0
    // 7 octobre 2004
    // Elmaestro
    //printf("\nlecture Newick : ");
    // TODO: Add your command handler code here
    //int FAIL =-1;
    int i,j,j1,k, a, a1, a2,a3 = 0, VertexNumber,numero = 0;
    char symbol, *string, *string1, *string2/* *string3,c ,**Et*/;
    int taxaPos; // le nombre de taxas recupéré
    int aretePos; // le nombre d'aretes recupéré
    char symbolOld =' ';
    int zz, xx,jj,ii;
    double longueur;
    char * tempString;
    int cpt=0;
    int na;
    char *string4 = (char*) malloc(1000);
    int temoin=0;

    a=0;

    //Correctness of the Newick format verification
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
        //printf("Unable to read your input data, please check it and try again...");
        //exit(-1);
    }

    j=0;
    do{
        symbol=newick[cpt++];
        if (symbol=='(') j++;
    }while(symbol != '\0');

    cpt=0;
    j1=0;

    do{
        symbol=newick[cpt++];
        if (symbol==')') j1++;
    }while(symbol != '\0');

    cpt=0;

    // verification des arités de l'arbre
    if (j1!=j) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); return -1;}
    //else if (j!=n-2) { printf("Incorrect Newick file format. Only trees with vertices of degree 1 and 3 are allowed by T-REX."); fclose (data); return -1;}

    k=0;

    do{
        symbol=newick[cpt++];
        if (symbol==',') k++;
    }while(symbol != '\0');

    cpt=0;
    //if (k!=(n-1)) { printf("Incorrect Newick file format. Number of objects must be equal to number of commas plus 1."); fclose (data); return -1;}

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

    if ((string == NULL)||(string1 == NULL)||(string2 == NULL)/*||(string3 == NULL)*/)
    { printf("Input data are too large or not a correct Newick file chosen"); return -1;}
    //printf("TEST1");
    a=0;

    do{
        symbol=newick[cpt++];
        if ((symbol!=' ')&&(symbol!='\n')&&(symbol!='\t')) { string[a++]=symbol; }
    }while(symbol !='\0');

    int boot_value;
    int temoin3 = 0;
    k=0; VertexNumber=n;
    //a1 = 0;
    //a2 = 0;
    taxaPos =1;    // nous allons commencer à mettre les taxas à la position 1
    aretePos = 1;
    boot_value=0;

    while (string[0] == '(')   // traiter toute la chaine
    {
        a1 = 0;
        a2 = 0;
        while( string[a2] != ')')  // traiter la paire () la plus profonde
        {
            if(string[a2] == '(') a1 = a2;  // retrouver ;a parenthèse ouvrante
            a2++;
        }


        // a   => contient la longueur de la chaine
        // a1  => contient le debut d'un noeud à traiter
        // a2  => contient la fin d'un noeud à traiter
        // a3  => délimite le noeud et sa longueur

        zz = a1+1;
        VertexNumber++;  // augmenter le nombre de noeuds
        boot_value=0;
        for ( ii = a1+1; ii <= a2; ii++)
        {// decortiquer cette chaine

            if (string[ii] == ':')
            {
                xx = 0;
                a3 = ii+1;

                if( string[zz] == '%')
                { // cela veut dire que c'est un  noeud que l'on traite

                    for ( jj = zz+1; jj < ii; jj++)
                    {
                        if(string[jj] == '|')
                            break;
                        string1[xx++] = string[jj];
                    }
                    temoin3=1;
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
                    // on recupère le nom du taxa

                    for(jj = zz; jj < ii; jj++)
                    {
                        lesNoms[taxaPos-1][xx++] = string[jj];
                    }
                    numero = taxaPos;
                    lesNoms[taxaPos-1][xx] = '\0';  // mettre la fin de chaine
                    taxaPos++;  // augmenter le nombre de taxas pris
                }

            }
            else if(string[ii] == ','  || string[ii] == ')')
            {
                xx = 0;
                zz = ii +1;   // faire pointer sur le prochain noeud
                for ( jj = a3; jj < ii; jj++)
                {
                    string1[xx++] = string[jj];
                }
                string1[xx++] = '\0';
                longueur = atof(string1);
                ARETE[aretePos++] = VertexNumber;
                ARETE[aretePos++] = numero;
                LONGUEUR[(aretePos-1)/ 2] = longueur;

                if(boot_value == 1){
                    //printf("\n%d--%d : %lf (%s)",VertexNumber,numero,longueur,string4);
                    boot_value=0;
                }
            }

        }

        // fin for pour traiter noeud
        //transcrire la nouvelle chaine
        xx = 0;
        for ( jj = 0; jj < (int)a1; jj++)
        {string2[xx++] = string[jj];}

        // ecrire le vertex
        itoa_(VertexNumber,string1,10);
        string2[xx++] = '%';   // indiquer que c'est un noeud
        for( jj = 0; jj < (int) strlen(string1); jj++)
        {string2[xx++] = string1[jj];}

        int temoin=0;

        // transcrire la fin
        for( jj = a2+1; jj <= a; jj++)  // il faut voir si c'est  <= a ou c < a
        {
            if(string[jj] != ':' && temoin==0){
                string2[xx++] = '|';
            }
            temoin = 1;
            //if(temoin==1)
            string2[xx++] = string[jj];

        }

        // remplacer string par string2 en remplacant les pointeurs

        tempString = string;
        string = string2;
        string2 = tempString;
        tempString = 0;
        a = xx;  // mettre la longueur à jour

    } // fin du while pour traiter toute la string



    int root_existance = -1;

    for(jj=0;jj<n;jj++){
        if (strcmp(lesNoms[jj],"Root") == 0)
            root_existance = jj;
    }

    ARETE[aretePos++] = 0;
    ARETE[aretePos++] = 0;

    /*for(i=1;i<=2*(na);i++){
        ARETE[i-1] = ARETE[i];
    }*/

    i=0;
    int cpt_branches=0;
    do{
        i++;cpt_branches++;
        //    printf("\n%d : %d",cpt_branches,ARETE[i]);
        ARETE[i-1] = ARETE[i];
        i++;
        //    printf("--%d",ARETE[i]);
        ARETE[i-1] = ARETE[i];
    }while(ARETE[i] != 0);

    for(i=1;i<=na;i++){
        LONGUEUR[i-1] = LONGUEUR[i];
    }

    //== recherche des branches connexes a celle de la racine;
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


    //=== on teste si il y a un noeud de degre 2 et un noeud de degre 1

    int * tableau = (int*)malloc((2*n)*sizeof(int));
    int deg2=-1,deg1=-1;
    for(i=1;i<=2*n;i++)
        tableau[i-1] = 0;
    for(i=1;i<=na;i++){
        tableau[ARETE[2*i-1]]++;
        tableau[ARETE[2*i-2]]++;
    }
//    printf("\n");
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

        //== modification de la branche pos1
        if(ARETE[2*pos1-1] == deg2)
            ARETE[2*pos1-1] = (ARETE[2*pos2-1] == deg2)?ARETE[2*pos2-2]:ARETE[2*pos2-1];
        else
            ARETE[2*pos1-2] = (ARETE[2*pos2-1] == deg2)?ARETE[2*pos2-2]:ARETE[2*pos2-1];

        //== suppression de la branche pos2
        for(i=pos2;i<=na;i++){
            LONGUEUR[i-1] = LONGUEUR[i];
            ARETE[2*i-1] = ARETE[2*(i+1)-1];
            ARETE[2*i-2] = ARETE[2*(i+1)-2];
        }
        na--;

    }

    (*kt) = 2*n-3 - (*kt);


    for(i = n; i >= 1; i--){
        strcpy(lesNoms[i],lesNoms[i-1]);
    }
    return n;
}


/**
 * utils functions created by Arthur Debeaupte
 * */

void writeInputTreeToFile(const std::string& filename, InputTree& tree) {
    std::ofstream file(filename);

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