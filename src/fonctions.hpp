//
//  fonctions.hpp
//  k means phylogenetic trees clustering
//
//  Created by Benjamin ALBERTELLI on 13/06/2022.
//

#ifndef fonctions_hpp
#define fonctions_hpp

#include <stdio.h>
#include "structures.h"
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <limits>
#include <string>
#include <cstring>
#include "utils_tree.hpp"

using namespace std;

void initInputTree(struct InputTree *);
void allocMemmory(struct InputTree *, int);
void freeInputTree(struct InputTree *, int);
void freeReducedTree(struct InputTree *, int);
void computeCriteria(double **, double **, int, struct CRITERIA *, double *, long int *, double *, long int*);
int nbSpeciesNewick(string);
void newickToMatrix(string, struct InputTree *);
int readInputFile(string, string, const char *, struct InputTree *, struct InputTree *, char *);
void InitCriteria(struct CRITERIA *, int);
void FreeCriteria(struct CRITERIA *, int);
void CreateSubStructures(struct InputTree *, int, int);
void help();
bool file_exists(const char *);
int readInput(int, const char *, struct InputTree *);
int readParameters(struct Parameters *);

#endif /* fonctions_hpp */
