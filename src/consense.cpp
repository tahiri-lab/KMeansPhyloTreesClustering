#include "consense.hpp"

/**
 * Initialise les options du programme de consensus et les fixe à leurs
 * valeurs par défaut.
 *
 * Dans la version actuelle, toutes les options interactives sont désactivées :
 * le mode MRe (Majority Rule Extended) est sélectionné par défaut, sans
 * enracinement, et l'arbre est écrit dans le fichier de sortie.
 *
 * Variables globales modifiées : ibmpc, ansi, didreroot, firsttree, spp,
 * col, tree_pairing, strict, mr, mre, ml, mlfrac, noroot, numopts, outgrno,
 * outgropt, trout, prntsets, progress, treeprint.
 */
void getoptions()
{
  long loopcount, loopcount2;
  char ch;
  boolean done, done1;

  /* Réglages initiaux */
  ibmpc          = IBMCRT;
  ansi           = ANSICRT;
  didreroot      = false;
  firsttree      = true;
  spp            = 0 ;
  col            = 0 ;

  /* Requis pour que les fonctions de cons.c fonctionnent correctement */
  tree_pairing   = NO_PAIRING ;

  strict = false;
  mr = false;
  mre = true;
  ml = false;
  mlfrac = 0.5;
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  prntsets = true;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
#ifdef WIN32
    phyFillScreenColor();
#endif
    ch = 'Y';
    done = true;
    if (!done) {
      if ((noroot && (ch == 'O')) || strchr("CRT1234",ch) != NULL) {
        switch (ch) {

        case 'C':
          if (strict) {
            strict = false;
            mr = true;
          } else {
            if (ml) {
              ml = false;
              mre = true;
            } else {
              if (mre) {
                mre = false;
                strict = true;
              } else {
                if (mr) {
                  mr = false;
                  ml = true;
                }
              }
            }
          }
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt) {
            numopts++;
            loopcount2 = 0;
            do {
#ifdef WIN32
              phyFillScreenColor();
#endif
              fflush(stdout);
              scanf("%ld%*[^\n]", &outgrno);
              getchar();
              done1 = (outgrno >= 1);
            countup(&loopcount2, 10);
            } while (done1 != true);
          }
          break;

        case 'R':
          noroot = !noroot;
          break;

        case 'T':
          initterminal(&ibmpc, &ansi);
          break;

        case '1':
          prntsets = !prntsets;
          break;

        case '2':
          progress = !progress;
          break;

        case '3':
          treeprint = !treeprint;
          break;

        case '4':
          trout = !trout;
          break;

        }
      }
    }
    countup(&loopcount, 100);
  } while (!done);
  if (ml) {
    do {
      fflush(stdout);
      scanf("%lf%*[^\n]", &mlfrac);
      getchar();
    } while ((mlfrac < 0.5) || (mlfrac > 1.0));
  }
}  /* getoptions */


/**
 * Parcourt la liste circulaire des frères d'un nœud pour vérifier sa
 * cohérence.
 *
 * Pour les feuilles (nœud NULL), retourne immédiatement. Sinon, parcourt
 * la chaîne next dans la limite de 1000 itérations jusqu'à revenir au
 * nœud de départ.
 *
 * Paramètre :
 *   p - pointeur vers le nœud à inspecter
 */
void count_siblings(node **p)
{
  node *tmp_node;
  int i;

  if (!(*p)) {
    /* Feuille : pas de frères */
    return;
  } else {
    tmp_node = (*p)->next;
  }

  for (i = 0 ; i < 1000; i++) {
    if (tmp_node == (*p)) {
      /* Tous les frères ont été parcourus */
      break;
    } else if (tmp_node) {
      tmp_node = tmp_node->next;
    } else {
      return ;
    }
  }
} /* count_siblings */


/**
 * Écrit récursivement l'arbre enraciné en p dans le fichier outtree au
 * format Newick.
 *
 * Pour les feuilles, le nom est écrit en remplaçant les espaces par '_'.
 * Pour les nœuds internes, les sous-arbres sont séparés par des virgules
 * entre parenthèses. Les longueurs de branches correspondent à la fréquence
 * d'apparition de chaque nœud (deltav), sauf en mode strict.
 *
 * Paramètre :
 *   p - nœud courant (racine du sous-arbre à écrire)
 */
void treeout(node *p)
{
  long i, n = 0;
  char c;
  node *q;
  double x;

  count_siblings (&p);

  if (p->tip) {
    /* Feuille : écrire le nom en remplaçant les espaces par '_' */
    for (i = 1; i <= MAXNCH; i++) {
      if (p->nayme[i - 1] != '\0')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = p->nayme[i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    col += n;
  } else {
    /* Nœud interne : écrire les sous-arbres entre parenthèses */
    putc('(', outtree);
    col++;
    q = p->next;
    while (q != p) {
      treeout(q->back);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      col++;
      if (col > 60) {
        putc('\n', outtree);
        col = 0;
      }
    }
    putc(')', outtree);
    col++;
  }

  if (p->tip)
    x = ntrees;
  else
    x = (double)p->deltav;

  if (p == root) {
    /* Fin de l'arbre */
    fprintf(outtree, ";\n");
    return;
  }

  /* Écrire la longueur de branche selon la magnitude de x */
  else {
    if (!strict) {
      if (x >= 100.0) {
        fprintf(outtree, ":%5.1f", x);
        col += 4;
      } else if (x >= 10.0) {
          fprintf(outtree, ":%4.1f", x);
          col += 3;
        } else if (x >= 1.00) {
            fprintf(outtree, ":%4.2f", x);
            col += 3;
          }
    }
  }
}  /* treeout */


/**
 * Fonction principale du module de consensus phylogénétique.
 *
 * Lit un fichier d'arbres au format Newick (argv[1]), calcule l'arbre
 * consensus selon les options choisies (MRe par défaut) et écrit le
 * résultat dans le fichier "outtree".
 *
 * Étapes :
 *   1. Ouverture du fichier d'arbres en entrée
 *   2. Initialisation des options via getoptions()
 *   3. Comptage des arbres et des feuilles
 *   4. Lecture des groupes de bipartitions (read_groups)
 *   5. Calcul du consensus (consensus)
 *   6. Écriture de l'arbre résultat (treeout)
 *   7. Libération de la mémoire et fermeture des fichiers
 *
 * Paramètre :
 *   argv - tableau d'arguments (argv[1] = chemin du fichier d'arbres)
 *
 * Retourne "0" en cas de succès.
 */
string main_consensus(char *argv[])
{
  pattern_elm  ***pattern_array;
  long trees_in = 0;
  long i, j;
  long tip_count = 0;
  node *p, *q;
#ifdef MAC
  argc = 1;
  argv[0] = "Consense";
#endif
  init(0, argv);
  const char *in_tree = static_cast<const char*> (argv[1]);
  cout<<argv[1]<<endl;
  /* Ouverture en mode binaire : ftell() est cassé sous WIN32 avec les fins de ligne UNIX */
  openfile(&intree, in_tree, "input tree file", "rb", argv[0], intreename);

  getoptions();
  ntrees = 0.0;
  maxgrp = 32767;   /* taille initiale de la table de hachage des groupes */
  lasti  = -1;

  const char *out_tree = static_cast<const char*> ("outtree");
  if (trout)
    openfile(&outtree, out_tree, "output tree file", "w", argv[0], outtreename);

  trees_in = countsemic(&intree);
  countcomma(&intree,&tip_count);
  tip_count++; /* countcomma donne le nombre brut de virgules, les feuilles = virgules + 1 */

  /* Lecture des groupes de bipartitions et calcul de l'arbre consensus */
  read_groups(&pattern_array, trees_in, tip_count, intree);
  nodep      = (pointarray)Malloc(2*(1+spp)*sizeof(node *));
  for (i = 0; i < spp; i++) {
    nodep[i] = (node *)Malloc(sizeof(node));
    for (j = 0; j < MAXNCH; j++)
      nodep[i]->nayme[j] = '\0';
    strncpy(nodep[i]->nayme, nayme[i], MAXNCH);
  }
  for (i = spp; i < 2*(1+spp); i++)
    nodep[i] = NULL;
  consensus(pattern_array, trees_in);
  printf("\n");
  if (trout) {
    treeout(root);
  }
  for (i = 0; i < spp; i++)
    free(nodep[i]);
  for (i = spp; i < 2*(1 + spp); i++) {
    if (nodep[i] != NULL) {
      p = nodep[i]->next;
      do {
        q = p->next;
        free(p);
        p = q;
      } while (p != nodep[i]);
      free(p);
    }
  }
  free(nodep);
  FClose(outtree);
  FClose(intree);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Consensus Program Done.\n\n");

#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif

  return "0";
}  /* main_consensus */
