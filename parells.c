/*************************************************************************/
// Fitxer: parells.c
// Autor: Xavier Sol Acha
// Descripci: Aqu estan les funcions per alinear parells de seqncies
/*************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "align.h"
#include "ent_seq.h"
#include "const.h"
#include "auxiliar.h"
#include "parells.h"
#include "malloc.h"

/* similarity *************************************************************/
// It calculates the similarity between the sequencies "seq1" and "seq2". 
// In the matrix "matrix" there are the similarities between every prefix of the two sequences.
float similitud(float **a, char *seq1, char *seq2, long longseq1, long longseq2)
{
    // use the global variables matpenal,alphabet,numsimb
    // SEQ1 and SEQ2 go from 0 to longseq-1

    float sum1, sum2, sum3, sum;
    int i, j;
    char c;

    a[0][0] = 0;

    for (i = 1; i <= longseq1; i++)
    {
        a[i][0] = a[i - 1][0] + matpenal[indexsimbol(seq1[i - 1])][numsimb];
        matriu_cami[i][0] = 'a';
    }

    for (j = 1; j <= longseq2; j++)
    {
        a[0][j] = a[0][j - 1] + matpenal[indexsimbol(seq2[j - 1])][numsimb];
        matriu_cami[0][j] = 'e';
    }

    for (i = 1; i <= longseq1; i++)
    {
        for (j = 1; j <= longseq2; j++)
        {
            // ve de dalt
            sum1 = a[i - 1][j] + matpenal[indexsimbol(seq1[i - 1])][numsimb];
            sum2 = a[i - 1][j - 1] + matpenal[indexsimbol(seq1[i - 1])][indexsimbol(seq2[j - 1])];
            sum3 = a[i][j - 1] + matpenal[indexsimbol(seq2[j - 1])][numsimb];
            a[i][j] = maxim_real(sum1, sum2, sum3, &c);
            matriu_cami[i][j] = c;
        }
    }
    return (a[longseq1][longseq2]);
}

/* alineament_optim
It returns the optimal alignment between the sequencies "seq1" and "seq2". 
It reserves the space in memory to make the alignment and the one that makes the call 
initial to recursive function.*/
