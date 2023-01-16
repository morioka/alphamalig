#include "align.h"
#include "const.h"
#include <stdio.h>
#include <stdlib.h>

void comprovarlecturaalfabet()
{
    int i, j;
    printf("\nNum simbols =%d\n", numsimb);

    for (i = 1; i <= numsimb; i++)
    printf("%c ", alfabet[i]);
    printf("\n");

    for (i = 1; i <= numsimb; i++)
    {
        for (j = 1; j <= i; j++)
        {
          printf("(%d,%d)=%f", i, j, matpenal[i][j]);
        }
        printf("\n");
    }
}

void comprovar_similitud(char *seq1, char *seq2, long longseq1, long longseq2)
{

    // char alin[2][2*MAXLONGSEQ];
    int i = longseq1;
    int j = longseq2;
    int pos = 0;
    char *alin0, *alin1;
    alin0 = (char *)malloc(2 * MAXLONGSEQ * sizeof(char));
    alin1 = (char *)malloc(2 * MAXLONGSEQ * sizeof(char));

    while ((i > 0) || (j > 0))
    {
        //    printf("cami[%d,%d]=%c,",i,j,matriu_cami[i][j]);
        if (matriu_cami[i][j] == 'a')
        {
            alin1[pos] = '-';
            alin0[pos] = seq1[i - 1];
            i--;
        }
        else if (matriu_cami[i][j] == 'd')
        {
            alin0[pos] = seq1[i - 1];
            alin1[pos] = seq2[j - 1];
            i--;
            j--;
        }
        else
        {
            alin1[pos] = seq2[j - 1];
            alin0[pos] = '-';
            j--;
        }
        //    printf("aparallem %c,%c", alin0[pos],alin1[pos]);
        pos++;
    }
    i = pos - 1;
    while (i > 80)
    {
        for (j = 0; j < 80; j++)
            printf("%c", alin0[i - j]);
        printf("\n");
        for (j = 0; j < 80; j++)
            printf("%c", alin1[i - j]);
        i = i - 80;
        printf("\n\n");
    }
    for (j = i; j >= 0; j--)
        printf("%c", alin0[j]);
    printf("\n");
    for (j = i; j >= 0; j--)
        printf("%c", alin1[j]);
    // for (i=pos-1;i>=0;i--) printf("%c",alin0[i]);
    printf("\n");
    // for(i=pos-1;i>=0;i--) printf("%c",alin1[i]);

    free(alin0);
    free(alin1);
}

void comprovar_carregar_sequencia(char *seq, long longseq)
{
    int i;
    printf("\n comprovacio carrega sequencia");
    for (i = 0; i < longseq; i++)
        printf("%c", seq[i]);
    printf("\n");
}

void comprovar_fitxer_temporal()
{
    char c;
    fitxer_temporal = fopen(nom_fitxer_temporal, "r");
    c = getc(fitxer_temporal);
    while (!feof(fitxer_temporal))
    {
        printf("%c", c);
        c = getc(fitxer_temporal);
    }
    fclose(fitxer_temporal);
}

void comprovar_matriu_similaritats(float **matriu)
{
    int i, j;
    printf("\n Similaritat entre les sequencies \n");
    for (i = 0; i < num_seqs - 1; i++)
    {
        for (j = i + 1; j < num_seqs; j++)
            printf("[%d,%d]=%f,", i + 1, j + 1, matriu[i][j]);
        printf("\n");
    }
}

void comprovar_info_cluster(float **info)
{
    int i, j;
    printf("\n Comprova info cluster\n");
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < numsimb; j++)
            printf("%f  ", info[i][j]);
        printf("\n");
    }
}

void comprova_matriu_cami(float **mat, int l1, int l2)
{
    int i, j;
    printf("l1=%d,l2=%d\n", l1, l2);
    for (i = 0; i <= l1; i++)
    {
        for (j = 0; j <= l2; j++)
            printf("%c,", matriu_cami[i][j]);
        for (j = 0; j <= l2; j++)
            printf("%f,", mat[i][j]);
        printf("\n");
    }
}
