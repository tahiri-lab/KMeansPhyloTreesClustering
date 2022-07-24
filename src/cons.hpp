#ifndef cons_hpp
#define cons_hpp

#define OVER              8
#define ADJACENT_PAIRS    1
#define CORR_IN_1_AND_2   2
#define ALL_IN_1_AND_2    3
#define NO_PAIRING        4
#define ALL_IN_FIRST      5
#define TREE1             8
#define TREE2             9

#define FULL_MATRIX       11
#define VERBOSE           22
#define SPARSE            33

/* Number of columns per block in a matrix output */
#define COLUMNS_PER_BLOCK 10

#include "phylip.hpp"

typedef struct pattern_elm {
  group_type *apattern;
  long *patternsize;
  double *length;
} pattern_elm;

/* begin hash table code */

#define NUM_BUCKETS 100

typedef struct namenode {
  struct namenode *next;
  plotstring naym;
  int hitCount;
} namenode;

typedef namenode **hashtype;

void read_groups (pattern_elm ****, long, long, FILE *);
void consensus(pattern_elm ***, long);

long namesGetBucket(plotstring);
void namesAdd(plotstring);
boolean namesSearch(plotstring);
void namesDelete(plotstring);
void namesClearTable(void);
void namesCheckTable(void);
void missingnameRecurs(node *p);

#ifndef OLDC
#define OLDC


/* prototypes */
void censor(void);
boolean compatible(long, long);
void elimboth(long);
void enterpartition (group_type*, long*);
void reorient(node* n);

/* function prototypes */
void initconsnode(node **, node **, node *, long, long, long *, long *,
                  initops, pointarray, pointarray, Char *, Char *, FILE *);
void   compress(long *);
void   sort(long);
void   eliminate(long *, long *);
void   printset(long);
void   bigsubset(group_type *, long);
void   recontraverse(node **, group_type *, long, long *);
void   reconstruct(long);
void   coordinates(node *, long *);
void   drawline(long i);

void   printree(void);
void   rehash(void);
void   enternodeset(node *r);
void   accumulate(node *);
void   dupname2(Char *, node *, node *);
void   dupname(node *);
void   missingname(node *);
void   gdispose(node *);
void   initreenode(node *);
void   reroot(node *, long *);

void store_pattern (pattern_elm ***, int);
boolean samename(naym, plotstring);
void reordertips(void);

void clean_up_final(void);
/* function prototypes */
#endif

extern long setsz;

#endif
