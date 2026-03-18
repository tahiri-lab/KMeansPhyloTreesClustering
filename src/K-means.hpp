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

extern bool withConsensus;

int main_kmeans(char **, vector <string>, double **, vector<int>, bool, int, int);

//Cleans up kmeans variables and memory allocations
void kmeans_cleanup(FILE *Output4, int kmax, int treeAmount, int **listr,
                    double *CHr, double *Wr,
                    int *list, int *no, int *howmany);

int setup_k_bounds(int k_max, int k_capacity, int k_min, int &k1, int &k2, int treeAmount, bool isBH, double *&CHr, double *&Wr);

//--Read the data
void ReadData1(int treeAmount,int &nmax,int &pmax);

//--Calculate the kmeans
void Assign(int &iran,int &treeAmount,int &nmax,int &k1,int* list,int* howmany,int* no,int &iassign,int random_number);

//--Squared distances to group centroids. Assign objects to nearest one
void Permute(int &treeAmount,int &nmax,int *iordre);
double f_RI(int Strouve[],int Sref[],int N);
double f_ARI(int Strouve[],int Sref[],const char *K_real,int group,int N);
void outStat(int Strouve[],int Sref[],char *criteria,int N,char *N_especes,char *percent,const char *K_real,int group,double score,/*int **listr,double *allScore,int k1, int k2,*/ vector <string> monTableau);

//SUPERTREES
double FO_super_tree(int &treeAmount,int &kmax,double** mat,int* list,int* howmany,double &SSE,int &currentK);

//fonctions sans passer par le centroid pour supertree
double DistanceCH(int &treeAmount,int &kmax,double** mat,int* list,double FO_new);

// CRITERE W
double FO_W(int &treeAmount,int &kmax,double** mat,int* list,int* howmany,double &SSE,int &currentK);
double DistanceW(int &treeAmount, int &kmax, int* list, double FO_new);

double arrondir(double num,int digits);
void conv2sameRef(int *Strouve,int *Sref, int N);

#endif /* K_means_hpp */

