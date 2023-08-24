//
//  utils_tree.hpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

#ifndef utils_tree_hpp
#define utils_tree_hpp

#include <stdio.h>
#include "structures.h"
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

double BipartitionDistance (int **, int ** ,int);
void Floyd(double **, double **, int, int);
void Floyd(double **, double **, double **, int, int);
void loadAdjacenceMatrix( double **, long int *, double *, int, int);
void odp1(double **, int *, int *, int *, int);
int Tree_edges (double **, long int *, double *, int, int);
int Bipartition_Table (double **, int **, int *, int);
int Table_Comparaison (int **, int **, int *, int *, int, int, int);
int findFils(double **, int, int);
struct TNoeud * CreerSousArbre(long int *, int *, double **, int, int);
void sortIntTab(int *, int, int);
void ParcoursArbre(struct TNoeud *, struct DescTree *);
void viderArbre(struct TNoeud *);
void AfficherArbre(struct TNoeud *, int);
static void xtoa (unsigned long, char *, unsigned, int);
char * itoa_(int, char *, int);
void filtrerMatrice(double **, double **, char **, char **, int, int, char *);
int ecrireMatrice(double **, const char *, int, char **);
void ajouterMatriceGene(double **, const char *, int, char **);
void TrierMatrices(double **, char **, char **, int);
int lectureNewick(string, long int *, double *, char **, int *);

/**
 * functions by Arthur Debeaupte
 * */
void writeInputTreeToFile(const std::string& filename, InputTree& tree);
void writeCriteriaToFile(const std::string& filename, CRITERIA& criteria);


#endif /* utils_tree_hpp */
