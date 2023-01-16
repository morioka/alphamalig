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

extern FILE *fitxer_entrada;  // Fitxer d'entrada de les seqncies
extern FILE *fitxer_temporal; // Fitxer temporal que emmagatzema les seqncies
extern FILE *fitxer_alfabet;  // Fitxer que emmagatzema l'alfabet en el formata adequat

extern char nom_fitxer_temporal[13]; // Nom del fitxer temporal

extern int numsimb;              // nombre de simbols de l'alfabet
extern char alfabet[MAXLONGALF]; // alfabet en posicions 1,2,..,numsimb
extern int num_seqs;             // Nombre de seqncies que hi ha al fitxer

extern float matpenal[MAXLONGALF][MAXLONGALF]; // PUNTUACIO match,mismatch i gap...
extern float **matriu_puntuacions;             // Matriu que contindr la puntuaci de l'alineament ptim entre totes les seqncies
extern char **matriu_cami;                     // Cont per on s'ha omplert cada posicio de la matriu ('e': esquerra, 'a':amunt, 'd': diagonal)

extern int args(int argc, char **argv);
extern void leer_alfabeto(FILE *fd);

#endif
