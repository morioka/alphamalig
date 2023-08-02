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

extern FILE *input_file;  // file of entry of the sequencies
extern FILE *temp_file; // temporary file that stores the sequencies
extern FILE *alphabet_file;  // File that stores the alphabet in the appropriate format

extern char temp_file_name[13]; // Temporary file name

extern int num_symbols;              // number of alphabet symbols
extern unsigned char alphabet[MAXLENALPHABET]; // alphabet in positions 1,2,..,num_symbols
extern int num_seqs;             // Number of sequences in the file

extern float matpenal[MAXLENALPHABET][MAXLENALPHABET]; // SCORE match,mismatch and gap...
extern float **score_matrix;             // Matrix containing the score of the ptim alignment between all the sequences
extern char **path_matrix;               // Count where each position of the matrix has been filled ('e': left, 'a':up, 'd': diagonal)

extern int args(int argc, char **argv);
extern void read_alphabet(FILE *fd);

#endif
