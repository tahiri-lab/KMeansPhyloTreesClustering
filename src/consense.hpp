#ifndef consense_hpp
#define consense_hpp

#include <stdio.h>

#include "cons.hpp"
#include <iostream>

using namespace std;

/* version 3.6. (c) Copyright 1993-2008 by the University of Washington.
   Written by Joseph Felsenstein, Hisashi Horino,
   Akiko Fuseki, Dan Fineman, Sean Lamont, and Andrew Keeffe.
   Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */


/* The following extern's refer to things declared in cons.c */

extern int tree_pairing;

extern char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
extern node *root;

extern long numopts, outgrno, col;
extern long maxgrp;               /* max. no. of groups in all trees found  */

extern boolean trout, firsttree, noroot, outgropt, didreroot, prntsets,
          progress, treeprint, goteof, strict, mr, mre, ml;
extern pointarray nodep;                 /* pointers to all nodes in tree */
extern group_type **grouping, **grping2, **group2;/* to store groups found  */
extern long **order, **order2, lasti;
extern group_type *fullset;
extern long tipy;

extern double trweight, ntrees, mlfrac;


void getoptions(void);
void count_siblings(node);
void treeout(node);
string main_consense(char *);

#endif /* consense_hpp */
