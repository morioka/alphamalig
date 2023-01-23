/*************************************************************************/
// File: align.c
// Author: Xavier Sol Acha
// Description: This file will contain the main program of the application
/*************************************************************************/

#include "align.h"
#include "auxiliar.h"
#include "ent_seq.h"
#include "parells.h"
#include "const.h"
#include "multiple.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <ctype.h>
#include "comprobacions.h"

FILE *input_file;  // file of entry of the sequencies
FILE *temp_file; // temporary file that stores the sequencies
FILE *alphabet_file;  // File that stores the alphabet in the appropriate format

char temp_file_name[13]; // Temporary file name

int num_symbols;              // number of alphabet symbols
char alphabet[MAXLONGALF]; // alphabet in positions 1,2,..,num_symbols
int num_seqs;             // Number of sequences in the file

float matpenal[MAXLONGALF][MAXLONGALF]; // SCORE match,mismatch and gap...
float **matriu_puntuacions;             // Matrix containing the score of the ptim alignment between all the sequences
char **matriu_cami;                     // Count where each position of the matrix has been filled ('e': left, 'a':up, 'd': diagonal)

int args(int argc, char **argv)
{
    // check the arguments
    int i;
    if (argc < 3)
    {
        fprintf(stderr, "Syntax: %s filealphabet filessequencesformatfasta \n", argv[0]);
        printf("The file format should be alphabet:\n");
        printf("    Number of symbols\n");
        printf("    Alphabet symbols separated by a white character\n");
        printf("    Gap penalty (there is only one type)\n");
        printf("and the following lines 1,2,...,i, contain 1 number, 2 numbers,...\n");
        printf("which are the scores of the symbol and the gap with the previous symbols");
        printf("EXAMPLE of 6 symbols (including the gap)\n");
        printf("6\n");
        printf("o p s n c -\n");
        printf(" 2     Score first symbol with same\n");
        printf("-1  2    Score second symbol with first and itself\n");
        printf("-2 -2  1   ...\n");
        printf("-2 -2  0  1\n");
        printf("-2 -2 -1 -1 1\n");
        printf("-2 -2  0  0 0 0\n");
        return (1);
    }
    else
        return (0);
}

// function that reads the alphabet, the score matrix and the gap
void read_alphabet(FILE *fd)
{
    int i, j;
    fscanf(fd, "%d\n", &num_symbols);

    for (i = 1; i <= num_symbols; i++)
    {
        if (num_symbols < 64)   // [A-Za-z0-9] + '-' (gap)
        {
            fscanf(fd, "%c ", &alphabet[i]);
        } else {
            int c;
            fscanf(fd, "%x ", &c);
            alphabet[i] = (char)c;
        }
    }
    //for (i = 1; i < num_symbols; i++)
    //{
    //    alphabet[i] = toupper(alphabet[i]);
    //}
    for (i = 1; i <= num_symbols; i++)
    {
        for (j = 1; j <= i; j++)
        {
            fscanf(fd, "%f", &matpenal[i][j]);
            matpenal[j][i] = matpenal[i][j];
        }
    }
}

int main(int argc, char *argv[])
{
    int i = 0, j, k, l, estat;
    long long_seq, long_seq2;
    char *seq, *seq2;
    float **matriu, res;
    char pref[6];
    int num_pref;

    if (args(argc, argv) == 1)
        exit(-1);
    if ((alphabet_file = fopen(argv[1], "r")) == NULL)
    {
        printf("Error reading the alphabet file\n");
        exit(1);
    }
    read_alphabet(alphabet_file);
    check_reading_alphabet();
    // initialize the global variables alphabet, matpenal, num_symbols
    if ((input_file = fopen(argv[2], "r")) == NULL)
    {
        printf("Error reading the sequence file\n");
        exit(1);
    }

    // We delete temporary files. We close the ones we don't need at the moment
    num_pref = generate_prefix_files(argv[1]);
    fill_string_prefix(pref, num_pref);
    sprintf(temp_file_name, "%s1.tmp", pref);
    sprintf(clusters_filename, "%s2.tmp", pref);
    sprintf(cluster_filename_1, "%s3.tmp", pref);
    sprintf(cluster_filename_2, "%s4.tmp", pref);
    temp_file = fopen(cluster_filename_1, "w");
    fclose(temp_file);
    temp_file = fopen(cluster_filename_2, "w");
    fclose(temp_file);
    temp_file = fopen(clusters_filename, "w");
    fclose(temp_file);
    temp_file = fopen(temp_file_name, "w");
    if (read_sequences_file() == -1)
    // have the sequences in the temporary file in a special format
    // number, name, length, sequence (numeroseq,nom,longitud,sequencia)
    {
        printf("Too many (>%d) or too few (<2) input sequences\n", MAXSEQ);
        remove(clusters_filename);
        remove(cluster_filename_1);
        remove(cluster_filename_2);
        remove(temp_file_name);
        exit(-1);
    }
    else // correct number of sequences
    {

        fclose(temp_file);
        //       check_temp_file();
        temp_file = fopen(temp_file_name, "r");

        estat = 0;
        matriu = (float **)malloc(MAXLONGALIN * sizeof(float *));
        if (matriu == NULL)
        {
          estat = -1;
          printf("Out of memory\n");
          remove(clusters_filename);
          remove(cluster_filename_1);
          remove(cluster_filename_2);
          remove(temp_file_name);
          exit(-1);
        }
        // there is enough memory for 2000 items
        i = 0;
        estat = 0;
        while ((i < MAXLONGALIN) && (estat == 0))
        {
          matriu[i] = (float *)malloc(MAXLONGALIN * sizeof(float));
          if (matriu[i] == NULL)
          {
            estat = -1;
            printf("Out of memory\n");
            remove(clusters_filename);
            remove(cluster_filename_1);
            remove(cluster_filename_2);
            remove(temp_file_name);
            exit(-1);
          }
          i++;
        }
        // there is enough memory for a 2000x2000 matrix

        estat = 0;
        matriu_cami = (char **)malloc(MAXLONGALIN * sizeof(char *));
        if (matriu_cami == NULL)
        {
          estat = -1;
          printf("Out of memory\n");
          remove(clusters_filename);
          remove(cluster_filename_1);
          remove(cluster_filename_2);
          remove(temp_file_name);
          exit(-1);
        }
        i = 0;
        estat = 0;
        while ((i < MAXLONGALIN) && (estat == 0))
        {
          matriu_cami[i] = (char *)malloc(MAXLONGALIN * sizeof(char));
          if (matriu_cami[i] == NULL)
          {
            estat = -1;
            printf("Out of memory\n");
            remove(clusters_filename);
            remove(cluster_filename_1);
            remove(cluster_filename_2);
            remove(temp_file_name);
            exit(-1);
          }
          i++;
        }

        matriu_puntuacions = (float **)malloc(num_seqs * sizeof(float *));
        if (matriu_puntuacions == NULL)
        {
          estat = -1;
        }
        i = 0;
        estat = 0;
        while ((i < num_seqs) && (estat == 0))
        {
            matriu_puntuacions[i] = (float *)malloc(num_seqs * sizeof(float));
            if (matriu_puntuacions[i] == NULL)
            {
              estat = -1;
              printf("Out of memory\n");
              remove(clusters_filename);
              remove(cluster_filename_1);
              remove(cluster_filename_2);
              remove(temp_file_name);
              exit(-1);
            }
            i++;
        }

        j = 0;
        seq = (char *)malloc(MAXLONGSEQ * sizeof(char));
        seq2 = (char *)malloc(MAXLONGSEQ * sizeof(char));

        // there is enough memory for all variables
        // printf("\n 1");
        while (j < (num_seqs - 1))
        {
            long_seq = dona_longitud_seq(j);
            carregar_sequencia(seq, j);
            // check_load_sequence(seq,long_seq);
            k = j + 1;
            printf("\n");
            while (k < num_seqs)
            {
                printf(".");
                long_seq2 = dona_longitud_seq(k);
                carregar_sequencia(seq2, k);
                // check_load_sequence(seq2,long_seq2);
                res = similarity(matriu, seq, seq2, long_seq, long_seq2);
                //      check_path_matrix(matriu,long_seq,long_seq2);
                //      check_similarity(seq, seq2, long_seq, long_seq2);
                matriu_puntuacions[j][k] = res;
                matriu_puntuacions[k][j] = res;
                k++;
            }
            j++;
        }
        check_similarity_matrix(matriu_puntuacions);
        free(seq);
        free(seq2);
        // printf("\n 2");
        multiple_alignment(matriu, atoi(argv[2]));
        // second argument is the output format
        // printf("\n 3");
        i = 0;
        while (i < num_seqs)
        {
            free(matriu_puntuacions[i]);
            i++;
        }
        free(matriu_puntuacions);

        i = 0;
        while (i < MAXLONGALIN)
        {
            free(matriu[i]);
            i++;
        }
        free(matriu);

        i = 0;
        while (i < MAXLONGALIN)
        {
            free(matriu_cami[i]);
            i++;
        }
        free(matriu_cami);
        fclose(temp_file);
        remove(clusters_filename);
        remove(temp_file_name);
        remove(cluster_filename_1);
        remove(cluster_filename_2);
    }
}
