#include "align.h"
#include "const.h"
#include <stdio.h>
#include <stdlib.h>

void check_reading_alphabet()
{
    int i, j;
    printf("\nNum of symbols =%d\n", num_symbols);

    for (i = 1; i <= num_symbols; i++)
        if ((alphabet[i] > 0x20) && (alphabet[i] < 0xff))   // printable?
        {
            printf("%c ", alphabet[i]);
        } else {
            printf("%02x ", alphabet[i]);
        }
    printf("\n");

    for (i = 1; i <= num_symbols; i++)
    {
        for (j = 1; j <= i; j++)
        {
          printf("(%d,%d)=%f", i, j, matpenal[i][j]);
        }
        printf("\n");
    }
}

void check_similarity(char *seq1, char *seq2, long longseq1, long longseq2)
{

    int i = longseq1;
    int j = longseq2;
    int pos = 0;
    char *alin0, *alin1;
    alin0 = (char *)malloc(2 * MAXLONGSEQ * sizeof(char));
    alin1 = (char *)malloc(2 * MAXLONGSEQ * sizeof(char));

    while ((i > 0) || (j > 0))
    {
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
    printf("\n");

    free(alin0);
    free(alin1);
}

void check_load_sequence(char *seq, long longseq)
{
    int i;
    printf("\n check load sequence");
    for (i = 0; i < longseq; i++)
        printf("%c", seq[i]);
    printf("\n");
}

void check_temp_file()
{
    char c;
    temp_file = fopen(temp_file_name, "r");
    c = getc(temp_file);
    while (!feof(temp_file))
    {
        printf("%c", c);
        c = getc(temp_file);
    }
    fclose(temp_file);
}

void check_similarity_matrix(float **matriu)
{
    int i, j;
    printf("\n Similarity between the sequences \n");
    for (i = 0; i < num_seqs - 1; i++)
    {
        for (j = i + 1; j < num_seqs; j++)
            printf("[%d,%d]=%f,", i + 1, j + 1, matriu[i][j]);
        printf("\n");
    }
}

void check_cluster_info(float **info)
{
    int i, j;
    printf("\n Check cluster info\n");
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < num_symbols; j++)
            printf("%f  ", info[i][j]);
        printf("\n");
    }
}

void check_path_matrix(float **mat, int l1, int l2)
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
