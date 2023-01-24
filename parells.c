/*************************************************************************/
// File: parells.c
// Author: Xavier Sol Acha
// Description: Here are the functions to align pairs of sequences
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
float similarity(float **a, char *seq1, char *seq2, long lenseq1, long lenseq2)
{
    // use the global variables matpenal,alphabet,num_symbols
    // SEQ1 and SEQ2 go from 0 to lenseq-1

    float sum1, sum2, sum3, sum;
    int i, j;
    char c;

    a[0][0] = 0;

    for (i = 1; i <= lenseq1; i++)
    {
        a[i][0] = a[i - 1][0] + matpenal[symbol_index(seq1[i - 1])][num_symbols];
        path_matrix[i][0] = 'a';
    }

    for (j = 1; j <= lenseq2; j++)
    {
        a[0][j] = a[0][j - 1] + matpenal[symbol_index(seq2[j - 1])][num_symbols];
        path_matrix[0][j] = 'e';
    }

    for (i = 1; i <= lenseq1; i++)
    {
        for (j = 1; j <= lenseq2; j++)
        {
            // it comes from above
            sum1 = a[i - 1][j] + matpenal[symbol_index(seq1[i - 1])][num_symbols];
            sum2 = a[i - 1][j - 1] + matpenal[symbol_index(seq1[i - 1])][symbol_index(seq2[j - 1])];
            sum3 = a[i][j - 1] + matpenal[symbol_index(seq2[j - 1])][num_symbols];
            a[i][j] = real_max(sum1, sum2, sum3, &c);
            path_matrix[i][j] = c;
        }
    }
    return (a[lenseq1][lenseq2]);
}

/* alineament_optim
It returns the optimal alignment between the sequencies "seq1" and "seq2". 
It reserves the space in memory to make the alignment and the one that makes the call 
initial to recursive function.*/
