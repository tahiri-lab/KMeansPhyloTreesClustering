//
//  K-means.hpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 14/06/2022.
//

#ifndef K_means_hpp
#define K_means_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>


using namespace std;

int main_kmeans(char **,vector <string>, double **, double **, double **, vector<int>, int, int *, int, int);

//--Read the data
void ReadData1(int &n,int &nmax,int &p,int &pmax,double** mat/* ,double* coord */,int* ishort,double* weight,double* colsum,int &ntran,char* nameb,int N);

//--Calculate the kmeans
void Assign(int &iran,int &n,int &nmax,int &k1,int &kmax,int* list,int* howmany,int* no,int &idebug,int &iassign,int &iseed, int random_number);

//--Squared distances to group centroids. Assign objects to nearest one
void CompSST(int &n,int &nmax,int &p,int &pmax,double** mat,double* weight,int* ishort,double &SST);
void Permute(int &iseed,int &n,int &nmax,int *iordre);
double f_RI(int Strouve[],int Sref[],int N);
double f_ARI(int Strouve[],int Sref[],const char *K_real,int group,int N);
void outStat(int Strouve[],int Sref[],char *criteria,int N,char *N_especes,char *percent,const char *K_real,int group,double score,int **listr,double *allScore,int k1, int k2 ,vector <string> monTableau);

//SUPERTREES
double FO_super_tree(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau);

//fonctions sans passer par le centroid pour supertree
double DistanceCH(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new, double facteur);

// CRITERE W
double FO_W(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau);
double DistanceW(int &n,int &kmax,double** mat, int* list, double** Ww, double FO_new, double facteur);

double arrondir(double num,int digits);
void conv2sameRef(int *Strouve,int *Sref, int N);

#endif /* K_means_hpp */

