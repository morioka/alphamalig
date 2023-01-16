/*************************************************************************/
// Fitxer: parells.c
// Autor: Xavier Sol� Acha
// Descripci�: Aqu� estan les funcions per alinear parells de seq��ncies
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

/* similitud *************************************************************/
// Calcula la similitud entre les seq��ncies "seq1" i "seq2". A la matriu
//"matriu" hi ha les similituds entre tot prefix de les dues seq��ncies.
float similitud(float **a, char *seq1, char *seq2, long longseq1, long longseq2)
{
    // usa les variables globals matpenal,alfabet,numsimb
    // seq1 i seq2 van de 0 a longseq-1

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
Retorna l'alineament �ptim entre les seq��ncies "seq1" i "seq2". �s la
que reserva l'espai a mem�ria per fer lalineament i la que fa la crida
inicial a la funci� recursiva.*/
