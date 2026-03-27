#include "cons.hpp"

/* =========================================================================
 * cons.cpp — Calcul du consensus d'arbres phylogénétiques
 *
 * Ce fichier implémente les algorithmes de construction d'arbres consensus
 * à partir d'un ensemble d'arbres phylogénétiques (format Newick).
 *
 * Les modes de consensus supportés sont :
 *   - strict : le clade doit être présent dans 100 % des arbres
 *   - mr     : règle de la majorité (> 50 % des arbres)
 *   - mre    : règle de la majorité élargie (tous les clades compatibles)
 *   - ml     : majorité avec seuil de fréquence personnalisé (mlfrac)
 * ========================================================================= */

int tree_pairing;

Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
node *root;

long numopts, outgrno, col, setsz;
long maxgrp;               /* nombre maximal de groupes distincts identifiés dans l'ensemble des arbres  */

/* Indicateurs du mode de consensus :
 *   mr     = règle de la majorité (le groupe doit apparaître dans plus de 50 % des arbres)
 *   mre    = règle de la majorité élargie (tous les groupes compatibles avec la majorité)
 *   ml     = règle de la majorité avec seuil de fréquence personnalisé mlfrac
 *   strict = consensus strict (le groupe doit apparaître dans 100 % des arbres) */
boolean trout, firsttree, noroot, outgropt, didreroot, prntsets,
               progress, treeprint, goteof, strict, mr=false, mre=false,
               ml=false; /* initialized all false for Treedist */
pointarray nodep;
pointarray treenode;

/* grouping : représentation sous forme de vecteur binaire de chaque clade identifié dans l'ensemble des arbres.
 * grping2, group2 : tableaux temporaires utilisés lors du rehashing et de l'élimination. */
group_type **grouping, **grping2, **group2;
double *lengths, *lengths2;
long **order, **order2, lasti;
group_type *fullset;
node *grbg;
long tipy;

double **timesseen, **tmseen2, **times2 ;
double *timesseen_changes, *tchange2;
double trweight, ntrees, mlfrac;

hashtype hashp;

/* =========================================================================
 * FONCTIONS DE TABLE DE HACHAGE POUR LA GESTION DES NOMS D'ESPECES
 * ========================================================================= */

/**
 * Calcule le seau de hachage associé à un nom.
 *
 * La fonction additionne les codes des caractères du nom, puis applique
 * un modulo sur le nombre de seaux disponibles.
 */
long namesGetBucket(plotstring searchname) {
  long i;
  long sum = 0;

  for (i = 0; (i < MAXNCH) && (searchname[i] != '\0'); i++) {
    sum += searchname[i];
  }
  return (sum % NUM_BUCKETS);
}

/**
 * Ajoute un nom dans la table de hachage.
 *
 * Le nom est inséré en tête de la liste chaînée du seau correspondant.
 * Aucun test de doublon n'est fait ici : si nécessaire, il faut appeler
 * namesSearch avant.
 */
void namesAdd(plotstring addname) {
  long bucket = namesGetBucket(addname);
  namenode *hp, *temp;

  temp = hashp[bucket];
  hashp[bucket] = (namenode *)Malloc(sizeof(namenode));

  hp = hashp[bucket];
  strcpy(hp->naym, addname);
  hp->next = temp;
  hp->hitCount = 0;
}

/**
 * Recherche un nom dans la table de hachage.
 *
 * Retourne true si le nom est trouvé, sinon false.
 * Si le nom est trouvé, son compteur hitCount est incrémenté.
 */
boolean namesSearch(plotstring searchname) {
  long i = namesGetBucket(searchname);
  namenode *p;

  p = hashp[i];
  if (p == NULL) {
    return false;
  }
  do {
    if (strcmp(searchname,p->naym) == 0) {
      p->hitCount++;
      return true;
    }
    p = p->next;
  } while (p != NULL);

  return false;
}

/**
 * Vérifie l'état de la table de hachage après les recherches de noms.
 *
 * - hitCount > 1 : le nom apparaît en double
 * - hitCount == 0 : le nom attendu n'a pas été retrouvé
 *
 * Après vérification, les compteurs sont remis à 0 pour une future réutilisation.
 */

void namesCheckTable(void) {
  namenode *p;
  long i;

  for (i=0; i< NUM_BUCKETS; i++) {
    p = hashp[i];
    while (p != NULL){
      if(p->hitCount >1){
        printf("\n\nERROR in user tree: duplicate name found: ");
        puts(p->naym);
        printf("\n\n");
        exxit(-1);
      } else if(p->hitCount == 0){
        printf("\n\nERROR in user tree: name %s not found\n\n\n",
               p->naym);
        exxit(-1);
      }
      p->hitCount = 0;
      p = p->next;
    }
  }
}

/**
 * Vide complètement la table de hachage.
 *
 * Chaque liste chaînée est libérée, puis le pointeur du seau est remis à NULL.
 */
void namesClearTable(void) {
  long i;
  namenode *p, *temp;

  for (i=0; i< NUM_BUCKETS; i++) {
    p = hashp[i];
    if (p != NULL) {
      do {
        temp = p;
        p = p->next;
        free(temp);
      } while (p != NULL);
    hashp[i] = NULL;
    }
  }
}
/* end hash table code */

/* =========================================================================
 * Node initialisation
 * ========================================================================= */

/**
 * initconsnode - initialise un nœud d'arbre lors de l'analyse Newick.
 *
 * Appelé par l'analyseur générique treeread() avec un sélecteur whichinit qui
 * indique quel aspect du nœud doit être configuré :
 *   bottom     - alloue et enregistre un nœud intérieur
 *   nonbottom  - alloue un nœud d'anneau non primaire
 *   tip        - alloue un nœud feuille et copie son nom
 *   length     - analyse et stocke la longueur de la branche
 *   treewt     - lit un poids d'arbre facultatif à partir d'un commentaire entre crochets
 *   unittrwt   - définit le poids de l'arbre à 1,0 (valeur par défaut lorsqu'aucun poids n'est fourni)
 *   hsnolength - marque la longueur de la branche comme manquante (-1)
 */
void initconsnode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str,
                        Char *ch, FILE *intree)
{
  /* initializes a node */
  long i;
  char c;
  boolean minusread;
  double valyew, divisor, fracchange;

  switch (whichinit) {
  case bottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->tip = false;
    for (i=0; i<MAXNCH; i++)
      (*p)->nayme[i] = '\0';
    nodep[(*p)->index - 1] = (*p);
    (*p)->v = 0;
    break;

  case nonbottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->v = 0;
    break;

  case tip:
    (*ntips)++;
    gnu(grbg, p);
    nodep[(*ntips) - 1] = *p;
    setupnode(*p, *ntips);
    (*p)->tip = true;
    strncpy ((*p)->nayme, str, MAXNCH);
    (*p)->v = 0;
    break;

  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    fracchange = 1.0;
    (*p)->v = valyew / divisor / fracchange;
    break;

  case treewt:
    if (!eoln(intree)) {
      if (fscanf(intree, "%lf", &trweight) == 1) {
        getch(ch, parens, intree);
        if (*ch != ']') {
          printf("\n\nERROR: Missing right square bracket\n\n");
          exxit(-1);
        } else {
          getch(ch, parens, intree);
          if (*ch != ';') {
            printf("\n\nERROR: Missing semicolon after square brackets\n\n");
            exxit(-1);
          }
        }
      }
      else {
        printf("\n\nERROR: Expecting tree weight in last comment field\n\n");
        exxit(-1);
      }
    }
    break;

  case unittrwt:
    /* This comes not only when setting trweight but also at the end of
     * any tree. The following code saves the current position in a
     * file and reads to a new line. If there is a new line then we're
     * at the end of tree, otherwise warn the user. This function should
     * really leave the file alone, so once we're done with 'intree'
     * we seek the position back so that it doesn't look like we did
     * anything */
    trweight = 1.0 ;
    i = ftell (intree);
    c = ' ';
    while (c == ' ')  {
      if (eoff(intree)) {
        fseek(intree,i,SEEK_SET);
        return;
      }
      c = gettc(intree);
    }
    fseek(intree,i,SEEK_SET);
    if ( c != '\n' && c!= '\r')
      printf("WARNING: Tree weight set to 1.0\n");
    if ( c == '\r' )
      if ( (c == gettc(intree)) != '\n')
        ungetc(c, intree);
    break;

  case hsnolength:
    (*p)->v = -1;         /* signal value that a length is missing */
    break;

  default:                /* cases hslength, iter, hsnolength      */
    break;                /* should there be an error message here?*/
  }
} /* initconsnode */


/**
 * Supprime les groupes trop rares pour figurer dans l'arbre consensus.
 *
 * Selon le mode choisi, un groupe est conservé seulement s'il satisfait
 * le critère de fréquence correspondant :
 *   - mre    : tous les groupes compatibles sont conservés (aucune suppression ici)
 *   - mr     : le groupe doit apparaître dans plus de 50 % des arbres
 *   - ml     : le groupe doit dépasser la fraction mlfrac
 *   - strict : le groupe doit être présent dans tous les arbres
 *
 * Les groupes ne satisfaisant pas le critère sont libérés et mis à NULL.
 */
void censor(void)
{
  long i;

  i = 1;
  do {
    if (timesseen[i-1])
      if (!(mre || (mr && (2*(*timesseen[i-1]) > ntrees))
                || (ml && ((*timesseen[i-1]) > mlfrac*ntrees))
                || (strict && ((*timesseen[i-1]) == ntrees)))) {
        free(grouping[i - 1]);
        free(timesseen[i - 1]);
        grouping[i - 1] = NULL;
        timesseen[i - 1] = NULL;
    }
    i++;
  } while (i < maxgrp);
} /* censor */


/**
 * Compacte le tableau des groupes en déplaçant tous les groupes non nuls
 * vers le début du tableau.
 *
 * Après censor(), certaines cases peuvent être NULL. Cette fonction les élimine
 * en décalant les groupes valides vers les premières positions. La valeur *n
 * est mise à jour pour indiquer le nombre de groupes valides restants.
 */
void compress(long *n)
{
  long i, j;

  i = 1;
  j = 1;
  do {
    while (grouping[i - 1] != NULL)
      i++;
    if (j <= i)
      j = i + 1;
    while ((grouping[j - 1] == NULL) && (j < maxgrp))
      j++;
    if (j < maxgrp) {
      grouping[i - 1] = (group_type *)Malloc(setsz * sizeof(group_type));
      timesseen[i - 1] = (double *)Malloc(sizeof(double));
      memcpy(grouping[i - 1], grouping[j - 1], setsz * sizeof(group_type));
      *timesseen[i - 1] = *timesseen[j - 1];
      free(grouping[j - 1]);
      free(timesseen[j - 1]);
      grouping[j - 1] = NULL;
      timesseen[j - 1] = NULL;
    }
  } while (j != maxgrp);
  (*n) = i - 1;
}  /* compress */


/**
 * Trie les groupes par ordre décroissant de fréquence (timesseen),
 * en maintenant la correspondance entre grouping[] et timesseen[].
 *
 * Algorithme : tri Shell (complexité O(n log² n) en moyenne).
 * Les deux tableaux grouping[] et timesseen[] sont triés simultanément.
 */
void sort(long n)
{
  long gap, i, j;
  group_type *stemp;
  double rtemp;

  gap = n / 2;
  stemp = (group_type *)Malloc(setsz * sizeof(group_type));
  while (gap > 0) {
    for (i = gap + 1; i <= n; i++) {
      j = i - gap;
      while (j > 0) {
        if (*timesseen[j - 1] < *timesseen[j + gap - 1]) {
          memcpy(stemp, grouping[j - 1], setsz * sizeof(group_type));
          memcpy(grouping[j - 1], grouping[j + gap - 1], setsz * sizeof(group_type));
          memcpy(grouping[j + gap - 1], stemp, setsz * sizeof(group_type));
          rtemp = *timesseen[j - 1];
          *timesseen[j - 1] = *timesseen[j + gap - 1];
          *timesseen[j + gap - 1] = rtemp;
        }
        j -= gap;
      }
    }
    gap /= 2;
  }
  free(stemp);
}  /* sort */


/**
 * Vérifie si deux groupes i et j sont compatibles dans l'arbre de consensus.
 *
 * Deux groupes (clades) représentés par leurs vecteurs binaires sont compatibles
 * si l'une des quatre conditions suivantes est vraie :
 *   1. Leur intersection est vide (ensembles disjoints)
 *   2. grouping[i] ⊆ grouping[j]  (i est inclus dans j)
 *   3. grouping[j] ⊆ grouping[i]  (j est inclus dans i)
 *   4. Leur union couvre le fullset (condition spéciale pour les arbres non enracinés)
 *
 * Retourne true si les deux groupes sont compatibles, false sinon.
 */
boolean compatible(long i, long j)
{
  long k;

  /* Condition 1 : intersection vide ? (grouping[i] ∩ grouping[j] == ∅) */
  boolean intersect_vide = true;
  for (k = 0; k < setsz; k++)
    if ((grouping[i][k] & grouping[j][k]) != 0) {
      intersect_vide = false;
      break;
    }
  if (intersect_vide) return true;

  /* Condition 2 : grouping[i] ⊆ grouping[j] ? */
  boolean i_inclus_j = true;
  for (k = 0; k < setsz; k++)
    if ((grouping[i][k] & ~grouping[j][k]) != 0) {
      i_inclus_j = false;
      break;
    }
  if (i_inclus_j) return true;

  /* Condition 3 : grouping[j] ⊆ grouping[i] ? */
  boolean j_inclus_i = true;
  for (k = 0; k < setsz; k++)
    if ((grouping[j][k] & ~grouping[i][k]) != 0) {
      j_inclus_i = false;
      break;
    }
  if (j_inclus_i) return true;

  /* Condition 4 (arbres non enracinés seulement) :
   * fullset \ grouping[i] \ grouping[j] == ∅ (union couvre tout) */
  if (!noroot) return false;
  for (k = 0; k < setsz; k++)
    if ((fullset[k] & ~grouping[i][k] & ~grouping[j][k]) != 0)
      return false;
  return true;
} /* compatible */


/**
 * Élimine les groupes incompatibles avec les groupes qui les précèdent
 * dans le tableau trié par fréquence décroissante.
 *
 * Pour chaque groupe i, on vérifie sa compatibilité avec tous les groupes
 * j < i déjà acceptés. Si une incompatibilité est détectée, le groupe i
 * est déplacé dans le tableau group2/times2 (groupes rejetés) et retiré
 * du tableau principal.
 *
 * Les compteurs *n et *n2 indiquent respectivement le nombre de groupes
 * conservés et le nombre de groupes rejetés.
 */
void eliminate(long *n, long *n2)
{
  long i, j, k;
  boolean comp;

  for (i = 2; i <= (*n); i++) {
    comp = true;
    for (j = 0; comp && (j <= i-2); j++) {
      if ((timesseen[j] != NULL) && *timesseen[j] > 0) {
        comp = compatible(i-1,j);
        if (!comp) {
          (*n2)++;
          times2[(*n2)-1] = (double *)Malloc(sizeof(double));
          group2[(*n2)-1] = (group_type *)Malloc(setsz * sizeof(group_type));
          *times2[(*n2)-1] = *timesseen[i - 1];
          memcpy(group2[(*n2) - 1], grouping[i - 1], setsz * sizeof(group_type));
          *timesseen[i - 1] = 0.0;
          for (k = 0; k < setsz; k++)
            grouping[i - 1][k] = 0;
        }
      }
    }
    if (*timesseen[i - 1] == 0.0) {
      free(grouping[i - 1]);
      free(timesseen[i -  1]);
      timesseen[i - 1] = NULL;
      grouping[i - 1] = NULL;
    }
  }
}  /* eliminate */


/**
 * Affichage des ensembles d'espèces retenus dans l'arbre consensus.
 *
 * Fonctionnalité désactivée dans cette version du projet :
 * l'écriture dans outfile a été retirée car la sortie texte n'est pas utilisée.
 */
void printset(long n)
{
  (void)n;
}  /* printset */


/**
 * Trouve le plus grand sous-ensemble de st parmi les n groupements existants.
 *
 * Utilisé lors de la reconstruction de l'arbre consensus pour identifier,
 * parmi tous les groupes, celui qui est le plus grand sous-ensemble strict
 * de st. Le résultat est copié dans st à la sortie.
 *
 * Paramètres :
 *   st - ensemble de départ (modifié en sortie)
 *   n  - nombre de groupements à considérer
 */
void bigsubset(group_type *st, long n)
{
  long i, j;
  group_type *su;
  boolean max, same;

  su = (group_type *)Malloc(setsz * sizeof(group_type));
  for (i = 0; i < setsz; i++)
    su[i] = 0;
  for (i = 0; i < n; i++) {
    max = true;
    for (j = 0; j < setsz; j++)
      if ((grouping[i][j] & ~st[j]) != 0)
        max = false;
    if (max) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i][j] != st[j])
          same = false;
      max = !same;
    }
    if (max) {
      for (j = 0; j < setsz; j++)
        if ((su[j] & ~grouping[i][j]) != 0)
          max = false;
      if (max) {
        same = true;
        for (j = 0; j < setsz; j++)
          if (su[j] != grouping[i][j])
            same = false;
        max = !same;
      }
      if (max)
        memcpy(su, grouping[i], setsz * sizeof(group_type));
    }
  }
  memcpy(st, su, setsz * sizeof(group_type));
  free(su);
}  /* bigsubset */


/**
 * Parcourt récursivement l'arbre pour ajouter le prochain nœud à l'arbre consensus.
 *
 * À chaque appel, la fonction décompose l'ensemble st en sous-ensembles maximaux
 * issus des groupements retenus, puis construit les nœuds intérieurs et feuilles
 * correspondants. Les pointeurs sont chaînés en anneau pour former la structure
 * d'arbre multifurcant.
 *
 * Paramètres :
 *   p        - nœud à créer (résultat)
 *   st       - ensemble d'espèces à placer sous ce nœud
 *   n        - nombre de groupements disponibles
 *   nextnode - prochain indice de nœud intérieur disponible
 */
void recontraverse(node **p, group_type *st, long n, long *nextnode)
{
  long i, j = 0, k = 0, l = 0;

  boolean found, same = 0, zero, zero2;
  group_type *tempset, *st2;
  node *q, *r;

  for (i = 1; i <= spp; i++) {  /* compte les espèces dans l'ensemble */
    if (i == ((l+1)*SETBITS+1)) l++;
    if (((1L << (i - 1 - l*SETBITS)) & st[l]) != 0) {
      k++;               /* k  est le nombre d'espèces dans l'ensemble */
      j = i;             /* j  pointe sur la dernière espèce de l'ensemble */
    }
  }
  if (k == 1) {           /* si une seule espèce : créer un nœud feuille */
    *p = nodep[j - 1];
    (*p)->tip = true;
    (*p)->index = j;
    return;
  }
  gnu(&grbg, p);          /* sinon : créer un nœud intérieur */
  (*p)->tip = false;
  (*p)->index = *nextnode;
  nodep[*nextnode - 1] = *p;
  (*nextnode)++;
  (*p)->deltav = 0.0;
  for (i = 0; i < n; i++) { /* cherche le groupement correspondant à st */
    same = true;
    for (j = 0; j < setsz; j++)
      if (grouping[i][j] != st[j])
        same = false;
    if (same)
      (*p)->deltav = *timesseen[i];
  }
  tempset = (group_type *)Malloc(setsz * sizeof(group_type));
  memcpy(tempset, st, setsz * sizeof(group_type));
  q = *p;
  st2 = (group_type *)Malloc(setsz * sizeof(group_type));
  memcpy(st2, st, setsz * sizeof(group_type));

  zero = true;      /* copie de l'ensemble faite ; vérifier si vide */
  for (j = 0; j < setsz; j++)
    if (tempset[j] != 0)
      zero = false;
  if (!zero)
    bigsubset(tempset, n);    /* trouver le plus grand sous-ensemble */
  zero = zero2 = false;       /* tempset est maintenant ce sous-ensemble */
  while (!zero && !zero2) {
    zero = zero2 = true;
    for (j = 0; j < setsz; j++) {
      if (st2[j] != 0)
        zero = false;
      if (tempset[j] != 0)
        zero2 = false;
    }
    if (!zero && !zero2) {
      gnu(&grbg, &q->next);
      q->next->index = q->index;
      q = q->next;
      q->tip = false;
      r = *p;
      recontraverse(&q->back, tempset, n, nextnode); /* placer sur l'arbre */
      *p = r;
      q->back->back = q;
      for (j = 0; j < setsz; j++)
        st2[j] &= ~tempset[j];    /* retirer le sous-ensemble de l'ensemble */
      memcpy(tempset, st2, setsz * sizeof(group_type));
      found = false;
      i = 1;
      while (!found && i <= n) {
        if (grouping[i - 1] != 0) {
          same = true;
          for (j = 0; j < setsz; j++)
            if (grouping[i - 1][j] != tempset[j])
              same = false;
        }
        if ((grouping[i - 1] != 0) && same)
          found = true;
        else
          i++;
      }
      zero = true;
      for (j = 0; j < setsz; j++)
        if (tempset[j] != 0)
          zero = false;
      if (!zero && !found)
        bigsubset(tempset, n);
    }
  }
  q->next = *p;
  free(tempset);
  free(st2);
}  /* recontraverse */


/**
 * Reconstruit l'arbre consensus à partir des groupements retenus.
 *
 * Lance la traversée récursive recontraverse() en partant du fullset
 * (ensemble complet de toutes les espèces) afin de construire l'arbre
 * depuis la racine.
 *
 * Paramètre :
 *   n - nombre de groupements valides à utiliser
 */
void reconstruct(long n)
{
  long nextnode;
  group_type *s;

  nextnode = spp + 1;
  s = (group_type *)Malloc(setsz * sizeof(group_type));
  memcpy(s, fullset, setsz * sizeof(group_type));
  recontraverse(&root, s, n, &nextnode);
  free(s);
}  /* reconstruct */


/**
 * Calcule les coordonnées (x, y) de chaque nœud pour l'affichage de l'arbre.
 *
 * Les coordonnées sont calculées récursivement :
 *   - Pour une feuille : xcoord = 0, ycoord = position verticale courante
 *   - Pour un nœud intérieur : xcoord = max(x enfants) + OVER,
 *     ycoord = milieu entre le premier et le dernier enfant
 *
 * Le paramètre tipy est incrémenté à chaque feuille visitée.
 */
void coordinates(node *p, long *tipy)
{
  node *q, *first, *last;
  long maxx;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  q = p->next;
  maxx = 0;
  while (q != p) {
    coordinates(q->back, tipy);
    if (!q->back->tip) {
      if (q->back->xcoord > maxx)
        maxx = (long)(q->back->xcoord);
    }
    q = q->next;
  }
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = maxx + OVER;
  p->ycoord = (long)((first->ycoord + last->ycoord) / 2);
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */



/**
 * Tente d'insérer une partition dans la liste des partitions connues.
 *
 * Si la partition s1 est déjà présente, elle n'est pas ajoutée.
 * Si elle est nouvelle, elle est insérée à la première place libre et
 * son compteur timesseen est initialisé à 1.
 *
 * Paramètres :
 *   s1 - vecteur binaire représentant la partition à insérer
 *   n  - nombre courant de partitions (incrémenté si ajout)
 */
void enterpartition (group_type *s1, long *n)
{
  long i, j;
  boolean found;

/* this stuff all to be rewritten but left here so pieces can be used */
  found = false;
  for (i = 0; i < (*n); i++) {  /* cherche si la partition est déjà présente */
    found = true;
    for (j = 0; j < setsz; j++) {  /* vérifie les deux parties de la partition */
      found = found && (grouping[i][j] == s1[j]);
      found = found && (group2[i][j] == (fullset[j] & (~s1[j])));
    }
    if (found)
      break;
  }
  if (!found) {    /* si absente, l'ajouter à la première case libre */
    grouping[i] = (group_type *)Malloc(setsz * sizeof(group_type));
    timesseen[i] = (double *)Malloc(sizeof(double));
    group2[i] = (group_type *)Malloc(setsz * sizeof(group_type));
    for (j = 0; j < setsz; j++)
      grouping[i][j] = s1[j];
    *timesseen[i] = 1;
    (*n)++;
  }
} /* enterpartition */


/**
 * Élimine les paires de groupes mutuellement incompatibles (cas Adams).
 *
 * Contrairement à eliminate(), qui garde les groupes compatibles avec tous
 * les précédents, cette fonction supprime les deux membres d'une paire
 * incompatible. Elle est utilisée dans le contexte du consensus Adams.
 */
void elimboth(long n)
{
  long i, j;
  boolean comp;

  for (i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      comp = compatible(i,j);
      if (!comp) {
        *timesseen[i] = 0.0;
        *timesseen[j] = 0.0;
      }
    }
    if (*timesseen[i] == 0.0) {
      free(grouping[i]);
      free(timesseen[i]);
      timesseen[i] = NULL;
      grouping[i] = NULL;
    }
  }
  if (*timesseen[n-1] == 0.0) {
    free(grouping[n-1]);
    free(timesseen[n-1]);
    timesseen[n-1] = NULL;
    grouping[n-1] = NULL;
  }
}  /* elimboth */


/**
 * Construit l'arbre consensus à partir des groupements accumulés.
 *
 * Étapes :
 *   1. censor()    — supprime les groupes trop rares
 *   2. compress()  — compacte le tableau
 *   3. sort()      — trie par fréquence décroissante (sauf mode strict)
 *   4. eliminate() — retire les groupes incompatibles
 *   5. reconstruct() — construit l'arbre à partir des groupes retenus
 *   6. coordinates() — calcule les coordonnées pour l'affichage
 *
 * Paramètres :
 *   pattern_array - tableau de motifs de partitions (utilisé par treedist)
 *   trees_in      - nombre total d'arbres traités
 */
void  consensus(pattern_elm ***pattern_array, long trees_in)
{
  long i, n, n2, tipy;

  group2 = (group_type **)  Malloc(maxgrp*sizeof(group_type *));
  for (i = 0; i < maxgrp; i++)
    group2[i] = NULL;
  times2 = (double **)Malloc(maxgrp*sizeof(double *));
  for (i = 0; i < maxgrp; i++)
    times2[i] = NULL;
  n2 = 0;
  censor();                /* supprimer les groupes trop rares */
  compress(&n);            /* compacter vers le début du tableau */
  if (!strict) {           /* éliminer les incompatibles (sauf mode strict) */
    sort(n);
    eliminate(&n, &n2);
    compress(&n);
    }
  reconstruct(n);
  tipy = 1;
  coordinates(root, &tipy);
  if (prntsets) {
    printset(n);
    /* Restaurer les groupes exclus dans grouping/timesseen pour affichage */
    for (i = 0; i < n2; i++) {
      if (!grouping[i]) {
        grouping[i] = (group_type *)Malloc(setsz * sizeof(group_type));
        timesseen[i] = (double *)Malloc(sizeof(double));
        }
      memcpy(grouping[i], group2[i], setsz * sizeof(group_type));
      *timesseen[i] = *times2[i];
    }
    n = n2;
    if (n2 != 0)
      printset(n);
  }
  free(nayme);
  for (i = 0; i < maxgrp; i++)
    free(grouping[i]);
  free(grouping);
  for (i = 0; i < maxgrp; i++)
    free(order[i]);
  free(order);
  for (i = 0; i < maxgrp; i++)
    if (timesseen[i] != NULL)
      free(timesseen[i]);
  free(timesseen);
}  /* consensus */


/**
 * Agrandit la table de hachage des groupements en doublant sa capacité.
 *
 * Lorsque la table est saturée (enternodeset ne trouve plus de place libre),
 * cette fonction alloue une nouvelle table de taille 2*maxgrp, y redistribue
 * tous les groupements existants via le même schéma de hachage (multiplication
 * par le ratio d'or conjugué), puis libère l'ancienne table.
 *
 * Le ratio d'or conjugué ((√5 - 1) / 2) est utilisé comme multiplicateur
 * irrrationnel pour minimiser les collisions.
 */
void rehash()
{
  group_type *s;
  long i, j;
  double partie_frac, valeur_hachage, ratio_dor;
  boolean done;

  long old_maxgrp = maxgrp;
  long new_maxgrp = maxgrp * 2;

  tmseen2 =     (double **)Malloc(new_maxgrp*sizeof(double *));
  grping2 = (group_type **)Malloc(new_maxgrp*sizeof(group_type *));
  order2 =        (long **)Malloc(new_maxgrp*sizeof(long *));
  lengths2 =     (double *)Malloc(new_maxgrp*sizeof(double));
  tchange2 =     (double *)Malloc(new_maxgrp*sizeof(double));

  for (i = 0; i < new_maxgrp; i++)
  {
    tmseen2[i] = NULL;
    grping2[i] = NULL;
    order2[i] = NULL;
    lengths2[i] = 0.0;
    tchange2[i] = 0.0;
  }

  ratio_dor = (sqrt(5.0) - 1) / 2;   /* ratio d'or conjugué : multiplicateur irrrationnel */
  s = (group_type *)Malloc(setsz * sizeof(group_type));

  for (i = 0; i < old_maxgrp; i++) {
    long old_index = *order[i];
    long new_index = -1;
    memcpy(s, grouping[old_index], setsz * sizeof(group_type));

    /* Calcule la valeur de hachage par somme pondérée des éléments du vecteur */
    valeur_hachage = 0.0;
    for (j = 0; j < setsz; j++)
      valeur_hachage += s[j] /* pow(2, SETBITS*j)*/;

    /* Applique le multiplicateur irrrationnel et extrait la partie fractionnaire */
    partie_frac = valeur_hachage * ratio_dor;
    new_index = (long)(new_maxgrp * (partie_frac - floor(partie_frac)));

    /* Sondage linéaire en cas de collision */
    done = false;
    while (!done) {
      if (!grping2[new_index])
      {
        grping2[new_index] = (group_type *)Malloc(setsz * sizeof(group_type));
        memcpy(grping2[new_index], grouping[old_index], setsz * sizeof(group_type));

        order2[i] = (long *)Malloc(sizeof(long));
        *order2[i] = new_index;

        tmseen2[new_index] = (double *)Malloc(sizeof(double));
        *tmseen2[new_index] = *timesseen[old_index];

        lengths2[new_index] = lengths[old_index];

        tchange2[new_index] = timesseen_changes[old_index];

        free(grouping[old_index]);
        free(timesseen[old_index]);
        free(order[i]);

        grouping[old_index] = NULL;
        timesseen[old_index] = NULL;
        order[i] = NULL;

        done = true; /* emplacement trouvé avec succès */

      } else {
        new_index++;
        if (new_index >= new_maxgrp) new_index -= new_maxgrp;
      }
    }
  }

  free(lengths);
  free(timesseen);
  free(grouping);
  free(order);
  free(timesseen_changes);

  free(s);

  /* Remplacement des anciens tableaux par les nouveaux */
  timesseen = tmseen2;
  grouping = grping2;
  lengths = lengths2;
  order = order2;
  timesseen_changes = tchange2;

  maxgrp = new_maxgrp;

}  /* rehash */


/**
 * Insère l'ensemble d'espèces d'un nœud dans la table de hachage des groupements.
 *
 * Si l'ensemble est déjà présent, son compteur timesseen est incrémenté.
 * S'il est absent et qu'une case est libre, il est inséré avec timesseen = trweight.
 * Si la table est pleine, rehash() est appelé pour l'agrandir, puis l'insertion
 * est relancée récursivement.
 *
 * Les ensembles complets (fullset) sont ignorés car ils ne correspondent pas
 * à des clades informatifs.
 *
 * Algorithme de hachage : multiplication par le ratio d'or conjugué +
 * sondage linéaire en cas de collision.
 */
void enternodeset(node* r)
{
  long i, j, start;
  double valeur_hachage;
  boolean done, same;
  double times;
  group_type *s;

  s = r->nodeset;

  /* Ne pas insérer les ensembles complets (non informatifs) */
  same = true;
  for (i = 0; i < setsz; i++)
    if (s[i] != fullset[i])
      same = false;
  if (same)
    return;

  times = trweight;

  /* Calcul du code de hachage par multiplication par le ratio d'or conjugué */
  valeur_hachage = 0.0;
  double ratio_dor = ((sqrt(5.0) - 1.0) / 2.0);
  for (i = 0; i < setsz; i++)
    valeur_hachage += s[i] * ratio_dor;

  /* Utilise la partie fractionnaire pour déterminer l'indice de départ */
  i = (long)(maxgrp * (valeur_hachage - floor(valeur_hachage))) + 1;
  start = i;
  done = false;

  while (!done) {
    if (grouping[i - 1]) {        /* si la case est occupée */
      same = false;
      if (!(timesseen[i-1] == 0)) {
        same = true;
        for (j = 0; j < setsz; j++) {
          if (s[j] != grouping[i - 1][j])
            same = false;
        }
      }
    }
    if (grouping[i - 1] && same) {   /* groupe déjà présent : incrémenter timesseen */
      *timesseen[i - 1] += times;
      lengths[i - 1] = nodep[r->index - 1]->v;
      done = true;
    } else if (!grouping[i - 1]) {   /* case vide : insérer le nouveau groupe */
      grouping[i - 1] = (group_type *)Malloc(setsz * sizeof(group_type));
      lasti++;
      order[lasti] = (long *)Malloc(sizeof(long));
      timesseen[i - 1] = (double *)Malloc(sizeof(double));
      memcpy(grouping[i - 1], s, setsz * sizeof(group_type));
      *timesseen[i - 1] = times;
      *order[lasti] = i - 1;
      done = true;
      lengths[i - 1] = nodep[r->index -1]->v;
    } else {                          /* collision : sondage linéaire */
      i++;
      if (i > maxgrp) i -= maxgrp;
    }
    if (!done && i == start) {  /* table saturée : agrandir et réessayer */
      rehash();
      done = true;
      enternodeset(r);
    }
  }
}  /* enternodeset */


/**
 * Parcourt récursivement l'arbre pour accumuler les ensembles d'espèces (nodesets).
 *
 * Pour chaque nœud, calcule nodeset comme le OU binaire des nodesets de tous
 * ses descendants. Pour une feuille, un seul bit est activé (correspondant
 * à son indice). Les nœuds non binaires (autre que l'anneau de 2) et les feuilles
 * sont ensuite insérés dans la table de hachage via enternodeset().
 *
 * Cette fonction est appelée après la lecture de chaque arbre.
 */
void accumulate(node *r)
{
  node *q;
  long i;

  /* Initialise nodeset à zéro ; alloue la mémoire la première fois rencontrée */
  if (!r->nodeset)
  {
    r->nodeset = (group_type *)Malloc(setsz * sizeof(group_type));
  }
  for (i = 0; i < setsz; i++)
  {
    r->nodeset[i] = 0L;
  }

  if (r->tip) {
    /* Feuille : active le bit correspondant à l'indice de l'espèce */
    i = (r->index-1) / (long)SETBITS;
    r->nodeset[i] = 1L << (r->index - 1 - i*SETBITS);
  }
  else {
    /* Nœud intérieur : OU binaire des ensembles de tous les enfants */
    for (q = r->next; q != r; q = q->next) {
      accumulate(q->back);
      for (i = 0; i < setsz; i++)
        r->nodeset[i] |= q->back->nodeset[i];
    }
  }

  /* Enregistrer le nodeset dans la table de hachage (sauf nœuds anneau binaires) */
  if ((!r->tip && (r->next->next != r)) || r->tip)
    enternodeset(r);
}  /* accumulate */


/**
 * Vérifie récursivement l'unicité des noms de feuilles dans l'arbre.
 *
 * Parcourt l'arbre en partant de p. Pour chaque feuille rencontrée, cherche
 * son nom dans la table de hachage hashp :
 *   - s'il est déjà présent → erreur (doublon)
 *   - sinon → l'ajouter à la table
 *
 * À la fin, tous les noms des feuilles sont enregistrés dans hashp.
 */
void dupname(node *p)
{
  node *q;

  if (p->tip) {
    if (namesSearch(p->nayme)) {
      printf("\n\nERROR in user tree: duplicate name found: ");
      puts(p->nayme);
      printf("\n\n");
      exxit(-1);
    } else {
      namesAdd(p->nayme);
    }
  } else {
    q = p;
    while (p->next != q) {
      dupname(p->next->back);
      p = p->next;
    }
  }
}  /* dupname */


/**
 * Vérifie récursivement que chaque nom de feuille est présent dans le premier arbre.
 *
 * Pour chaque feuille, recherche son nom dans la table de hachage (remplie par
 * dupname() lors de la lecture du premier arbre). Si un nom est absent, une
 * erreur fatale est signalée.
 */
void missingnameRecurs(node *p)
{
  node *q;

  if (p->tip) {
    if (!namesSearch(p->nayme)) {
      printf("\n\nERROR in user tree: name %s not found in first tree\n\n\n",
             p->nayme);
      exxit(-1);
    }
  } else {
    q = p;
    while (p->next != q) {
      missingnameRecurs(p->next->back);
      p = p->next;
    }
  }
}  /* missingnameRecurs */

/**
 * Point d'entrée pour la vérification des noms manquants.
 *
 * Lance missingnameRecurs() pour vérifier que tous les noms de feuilles
 * existent bien dans le premier arbre, puis appelle namesCheckTable() pour
 * s'assurer qu'aucun nom n'a été oublié.
 */
void missingname(node *p){
  missingnameRecurs(p);
  namesCheckTable();
} /* missingname */


/**
 * Libère récursivement tous les nœuds d'un arbre en les remettant dans
 * la liste de recyclage grbg.
 *
 * Parcourt l'arbre en profondeur d'abord. Pour les feuilles, chuck() est
 * appelé directement. Pour les nœuds intérieurs, les nœuds d'anneau sont
 * libérés avant le nœud principal.
 */
void gdispose(node *p)
{
  node *q, *r;

  if (p->tip) {
    chuck(&grbg, p);
    return;
  }
  q = p->next;
  while (q != p) {
    gdispose(q->back);
    r = q;
    q = q->next;
    chuck(&grbg, r);
  }
  chuck(&grbg, p);
}  /* gdispose */


/**
 * Parcourt l'arbre et copie les noms d'espèces des feuilles dans le tableau nayme[].
 *
 * Pour chaque nœud feuille rencontré, son nom (p->nayme) est copié à la position
 * correspondant à son indice (p->index - 1) dans le tableau global nayme[].
 */
void initreenode(node *p)
{
  node *q;

  if (p->tip) {
    memcpy(nayme[p->index - 1], p->nayme, MAXNCH);
  } else {
    q = p->next;
    while (q && q != p) {
      initreenode(q->back);
      q = q->next;
    }
  }
} /* initreenode */


/**
 * Réoriente récursivement les nœuds intérieurs de l'arbre après un réenracinement.
 *
 * Pour chaque nœud intérieur, s'assure que nodep[index-1] pointe bien sur le
 * nœud primaire de l'anneau. Si ce n'est pas le cas, met à jour nodep et
 * copie la longueur de branche depuis le nœud opposé (back).
 */
void reorient(node* n) {
  node* p;

  if ( n->tip ) return;
  if ( nodep[n->index - 1] != n )  {
    nodep[n->index - 1] = n;
    if ( n->back )
      n->v = n->back->v;
  }

  for ( p = n->next ; p != n ; p = p->next)
    reorient(p->back);
}


/**
 * Réenracine l'arbre en plaçant la racine au niveau du groupe externe (outgroup).
 *
 * Si la racine a exactement 2 enfants, les deux branches sont fusionnées et
 * la racine est déplacée vers l'outgroup sans créer de nouveaux nœuds.
 * Si la racine a 3 enfants ou plus, deux nouveaux nœuds d'anneau sont créés
 * pour former la nouvelle racine bifurquante.
 *
 * Après repositionnement, reorient() est appelé pour mettre à jour nodep[].
 *
 * Paramètres :
 *   outgroup - nœud qui sera adjacent à la nouvelle racine
 *   nextnode - prochain indice de nœud intérieur disponible (mis à jour)
 */
void reroot(node *outgroup, long *nextnode)
{
  long i;
  node *p, *q;
  double newv;

  /* Compter les enfants de la racine actuelle et trouver le dernier */
  p = root;
  i = 0;
  while (p->next != root) {
    p = p->next;
    i++;
  }
  if (i == 2) {
    /* Racine à 2 enfants : fusion des branches */
    q = root->next;

    newv = q->back->v + p->back->v;

    /* Si l'outgroup est déjà à la racine : réorganiser l'ordre des enfants */
    if (outgroup == p->back) {
      root->next = p;
      p->next = q;
      q->next = root;

      q->back->v = newv;
      p->back->v = 0;
      return;
    }
    if (outgroup == q) {
      p->back->v = newv;
      q->back->v = 0;
      return;
    }

    /* Détacher la racine en reliant ses deux enfants directement */
    q->back->back = p->back;
    p->back->back = q->back;
    p->back->v = newv;
    q->back->v = newv;
  } else { /* 3 enfants ou plus */
    p->next = root->next;              /* joindre les anciens nœuds de la racine */
    nodep[root->index-1] = root->next; /* root->next devient le nœud primaire */

    /* Créer les deux nouveaux nœuds d'anneau pour la nouvelle racine */
    gnu(&grbg, &root->next);
    q = root->next;
    gnu(&grbg, &q->next);
    p = q->next;
    p->next = root;
    q->tip = false;
    p->tip = false;
    nodep[*nextnode] = root;
    (*nextnode)++;
    root->index = *nextnode;
    root->next->index = root->index;
    root->next->next->index = root->index;
  }
  newv = outgroup->v;
  /* root est composé de 3 nœuds "flottants" */
  /* q == root->next, p == root->next->next */

  /* Attacher la racine au niveau de l'outgroup */
  q->back = outgroup;
  p->back = outgroup->back;
  outgroup->back->back = p;
  outgroup->back = q;
  outgroup->v = 0;
  outgroup->back->v = 0;
  root->v = 0;
  p->v = newv;
  p->back->v = newv;
  reorient(root);
}  /* reroot */


/**
 * Stocke les groupements du dernier arbre lu dans le tableau pattern_array.
 *
 * Pour chaque groupement présent dans grouping[] et apparu pour la première
 * fois depuis le dernier appel (timesseen > timesseen_changes), copie le
 * vecteur binaire dans pattern_array et met à jour timesseen_changes.
 *
 * Le champ patternsize est mis à jour pour indiquer le nombre de groupements
 * enregistrés pour cet arbre.
 *
 * Note : utilisé par treedist pour calculer la distance de Robinson-Foulds.
 *
 * Paramètres :
 *   pattern_array   - tableau de motifs à remplir
 *   trees_in_file   - indice de l'arbre courant dans le fichier
 */
void store_pattern (pattern_elm ***pattern_array, int trees_in_file)
{
  long i, total_groups=0, j=0, k;

  /* Compter le nombre de groupements présents dans l'arbre courant */
  for (i = 0 ; i < maxgrp ; i++)
    if ((grouping[i] != NULL) &&
       (*timesseen[i] > timesseen_changes[i]))
      total_groups++ ;

  /* Allouer les structures pour stocker les motifs binaires */
  for (i = 0 ; i < setsz ; i++) {
    pattern_array[i][trees_in_file]
      = (pattern_elm *) Malloc(sizeof(pattern_elm)) ;
    pattern_array[i][trees_in_file]->apattern =
      (group_type *) Malloc (total_groups * sizeof (group_type)) ;
    pattern_array[i][trees_in_file]->length =
      (double *) Malloc (maxgrp * sizeof (double)) ;
      for ( j = 0 ; j < maxgrp ; j++ ) {
        pattern_array[i][trees_in_file]->length[j] = -1;
      }
    pattern_array[i][trees_in_file]->patternsize = (long *)Malloc(sizeof(long));
  }
  j = 0;

  /* Copier chaque groupement nouveau dans pattern_array */
  for (i = 0 ; i < maxgrp ; i++)
    if (grouping[i] != NULL) {
      if (*timesseen[i] > timesseen_changes[i]) {
        for (k = 0 ; k < setsz ; k++)
          pattern_array[k][trees_in_file]->apattern[j] = grouping[i][k] ;
        pattern_array[0][trees_in_file]->length[j] = lengths[i];
        j++ ;

        timesseen_changes[i] = *timesseen[i] ;
          /*
             EWFIX.BUG.756

             updates timesseen_changes to the current value
             pointed to by timesseen

             treedist uses this to determine if group i has been seen
             by comparing timesseen_changes[i] (the count now) with
             timesseen[i] (the count after reading next tree)

             We could make treedist more efficient by not keeping
             timesseen (and groupings, etc) around, but doing it
             this way allows us to share code between treedist and
             consense.

          */
      }
    }
  *pattern_array[0][trees_in_file]->patternsize = total_groups;

}  /* store_pattern */


/**
 * Compare deux noms d'espèces et retourne true s'ils sont identiques.
 *
 * Paramètres :
 *   name1 - premier nom (type naym, tableau de caractères)
 *   name2 - second nom (type plotstring)
 */
boolean samename(naym name1, plotstring name2)
{
  return !(strncmp(name1, name2, MAXNCH));
}  /* samename */


/**
 * Réordonne nodep[] et les indices des feuilles pour correspondre à l'ordre
 * des espèces établi lors de la lecture du premier arbre.
 *
 * Pour chaque position i dans nayme[], cherche dans nodep[] la feuille dont
 * le nom correspond, puis échange les pointeurs et remet les indices à jour.
 *
 * Pré-condition : nayme[] contient les noms dans l'ordre du premier arbre,
 * nodep[] contient les feuilles dans l'ordre de lecture de l'arbre courant.
 */
void reordertips()
{
  long i, j;
  node *t;

  for (i = 0; i < spp-1; i++) {
    for (j = i + 1; j < spp; j++) {
      if (samename(nayme[i], nodep[j]->nayme)) {
        /* Échanger les pointeurs dans nodep[] et mettre à jour les indices */
        t = nodep[i];

        nodep[i] = nodep[j];
        nodep[i]->index = i+1;

        nodep[j] = t;
        nodep[j]->index = j+1;

        break;  /* passer au prochain i */
      }
    }
  }
}  /* reordertips */


/**
 * Lit tous les arbres du fichier et accumule les groupements dans la table de hachage.
 *
 * Pour chaque arbre :
 *   1. Alloue les nœuds et lit l'arbre au format Newick (treeread)
 *   2. Pour le premier arbre : vérifie les doublons de noms (dupname) et
 *      initialise le tableau nayme[] (initreenode)
 *   3. Pour les arbres suivants : vérifie que tous les noms sont présents
 *      (missingname) et réordonne les feuilles (reordertips)
 *   4. Réenracine si nécessaire (reroot)
 *   5. Accumule les ensembles d'espèces dans la table de hachage (accumulate)
 *   6. Libère la mémoire de l'arbre (gdispose)
 *   7. Si tree_pairing actif : stocke le motif de l'arbre (store_pattern)
 *
 * Paramètres :
 *   pattern_array - tableau de motifs à remplir (modifié)
 *   total_trees   - nombre maximum d'arbres attendus
 *   tip_count     - nombre de feuilles par arbre
 *   intree        - fichier d'entrée contenant les arbres Newick
 */
void read_groups (pattern_elm ****pattern_array,
        long total_trees, long tip_count, FILE *intree)
{
  int i, j, k;
  boolean haslengths, initial;
  long nextnode, trees_read = 0;

  /* Allocation des structures principales *******************************/
  grouping  = (group_type **)  Malloc(maxgrp*sizeof(group_type *));
  lengths  = (double *)  Malloc(maxgrp*sizeof(double));

  timesseen_changes = (double*)Malloc(maxgrp*sizeof(double));
  for (i = 0; i < maxgrp; i++)
    timesseen_changes[i] = 0.0;

  for (i = 0; i < maxgrp; i++)
    grouping[i] = NULL;

  order     = (long **) Malloc(maxgrp*sizeof(long *));
  for (i = 0; i < maxgrp; i++)
    order[i] = NULL;

  timesseen = (double **)Malloc(maxgrp*sizeof(double *));
  for (i = 0; i < maxgrp; i++)
    timesseen[i] = NULL;

  nayme = (naym *)Malloc(tip_count*sizeof(naym));
  hashp = (hashtype)Malloc(sizeof(namenode) * NUM_BUCKETS);
  for (i=0;i<NUM_BUCKETS;i++) {
      hashp[i] = NULL;
  }

  /* setsz : nombre de mots nécessaires pour représenter tip_count espèces en bits */
  setsz = (long)ceil((double)tip_count/(double)SETBITS);

  if (tree_pairing != NO_PAIRING)
    {
      /* Allocation de pattern_array maintenant que setsz est connu */
      (*pattern_array) =
        (pattern_elm ***)Malloc(setsz * sizeof(pattern_elm **));

      for (j = 0 ; j < setsz ; j++)
      {
          (*pattern_array)[j] =
            (pattern_elm **)Malloc(total_trees * sizeof(pattern_elm *));
        for(k = 0 ; k < total_trees ; k++ )
        {
            (*pattern_array)[j][k] = NULL;
        }
      }
    }

  /* Construction du fullset : tous les bits des tip_count espèces activés */
  fullset = (group_type *)Malloc(setsz * sizeof(group_type));
  for (j = 0; j < setsz; j++)
    fullset[j] = 0L;
  k = 0;
  for (j = 1; j <= tip_count; j++) {
    if (j == ((k+1)*SETBITS+1)) k++;
    fullset[k] |= 1L << (j - k*SETBITS - 1);
  }
  /* Fin de l'allocation **************************************************/

  firsttree = true;
  grbg = NULL;
  initial = true;
  while (!eoff(intree)) {          /* lire jusqu'à la fin du fichier d'arbres */
    for (i = 0; i < maxgrp; i++) {
      lengths[i] = -1;
    }
    goteof = false;
    nextnode = 0;
    haslengths = true;
    allocate_nodep(&nodep, &intree, &spp);
    assert(spp == tip_count);
    treeread(intree, &root, treenode, &goteof, &firsttree, nodep,
              &nextnode, &haslengths, &grbg, initconsnode,true,-1);
    if (!initial) {
      missingname(root);
      reordertips();
    } else {
      initial = false;
      dupname(root);
      initreenode(root);
    }
    if (goteof)
      continue;
    ntrees += trweight;
    if (noroot) {
      reroot(nodep[outgrno - 1], &nextnode);
      didreroot = outgropt;
    }
    accumulate(root);
    gdispose(root);
    for (j = 0; j < 2*(1+spp); j++)
      nodep[j] = NULL;
    free(nodep);
    if (tree_pairing != NO_PAIRING) {
        /* Stocker le motif binaire de cet arbre pour le calcul des distances */
      store_pattern ((*pattern_array), trees_read) ;
      trees_read++ ;
    }
  }
} /* read_groups */


/**
 * Libère toute la mémoire allouée lors du traitement des arbres.
 *
 * Libère les tableaux grouping[], order[], timesseen[], ainsi que les
 * structures globales (nayme, timesseen_changes, fullset, lengths).
 * Vide et libère également la table de hachage des noms (hashp).
 *
 * À appeler en fin de traitement pour éviter les fuites mémoire.
 */
void clean_up_final()
{
    long i;
    for(i=0;i<maxgrp;i++)
    {
        if(grouping[i] != NULL) {
            free(grouping[i]);
        }
        if(order[i] != NULL) {
            free(order[i]);
        }
        if(timesseen[i] != NULL) {
            free(timesseen[i]);
        }
    }
    free(grouping);
    free(nayme);
    free(order);
    free(timesseen);
    free(timesseen_changes);
    free(fullset);
    free(lengths);

    namesClearTable();
    free(hashp);
}
