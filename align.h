/*************************************************************************/
// Fitxer: align.h
// Autor: Xavier Sol Acha
// Descripci: Aquest fitxer contindr totes les estructures de dades del
// programa, aix com les diverses variables globals que necessitem.
/*************************************************************************/
#include <stdio.h>
#include "const.h"

#ifndef _ALIGN_H_
#define _ALIGN_H_

extern FILE *fitxer_entrada;  // file of entry of the sequencies
extern FILE *fitxer_temporal; // temporary file that stores the sequencies
extern FILE *fitxer_alfabet;  // File that stores the alphabet in the appropriate format

extern char nom_fitxer_temporal[13]; // Temporary file name

extern int numsimb;              // number of alphabet symbols
extern char alfabet[MAXLONGALF]; // alphabet in positions 1,2,..,numsimb
extern int num_seqs;             // Number of sequences in the file

extern float matpenal[MAXLONGALF][MAXLONGALF]; // SCORE match,mismatch and gap...
extern float **matriu_puntuacions;             // Matrix containing the score of the ptim alignment between all the sequences
extern char **matriu_cami;                     // Count where each position of the matrix has been filled ('e': left, 'a':up, 'd': diagonal)

extern int args(int argc, char **argv);
extern void leer_alfabeto(FILE *fd);

#endif
