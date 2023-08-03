/*************************************************************************/
// File: align.c
// Author: Xavier Sol Acha
// Description: This file will contain the main program of the application
/*************************************************************************/

#include "align.h"
#include "auxiliary.h"
#include "ent_seq.h"
#include "pairs.h"
#include "const.h"
#include "multiple.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <ctype.h>
#include "checks.h"
#include <assert.h>
#include <string.h>

FILE *input_file;  // file of entry of the sequencies
FILE *temp_file; // temporary file that stores the sequencies
FILE *alphabet_file;  // File that stores the alphabet in the appropriate format

char temp_file_name[13]; // Temporary file name

int num_symbols;              // number of alphabet symbols
unsigned char alphabet[MAXLENALPHABET]; // alphabet in positions 1,2,..,num_symbols
int num_seqs;             // Number of sequences in the file

float matpenal[MAXLENALPHABET][MAXLENALPHABET]; // SCORE match,mismatch and gap...
float **score_matrix;             // Matrix containing the score of the ptim alignment between all the sequences
char **path_matrix;                     // Count where each position of the matrix has been filled ('e': left, 'a':up, 'd': diagonal)

int hex_output_mode;


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

void check_alphabet()
{
    int i;

    // prohibited character
    for (i = 1; i <= num_symbols; i++)
    {
        assert(alphabet[i] != 0x00);    // # NUL (0x00)
        assert(alphabet[i] != 0x3e);    // # > (0x3e)
        assert(alphabet[i] != 0x3d);    // # = (0x3d)
        assert(alphabet[i] != 0x3c);    // # < (0x3c)
        assert(alphabet[i] != 0x20);    // # Space (0x20)
        assert(alphabet[i] != 0x0d);    // # Carriage Return (0x0d)
        assert(alphabet[i] != 0x0a);    // # Line Feed (0x0a)
    }

    // must-have-at-tail character (Gap)
    assert(alphabet[num_symbols] == 0x2d);    // # - (0x2d)

    // check duplicates
    int check_alphabet[256];
    for (i = 0; i < 256; i++)
        check_alphabet[i] = 0;
    for (i = 1; i <= num_symbols; i++)
        check_alphabet[alphabet[i]] += 1;
    for (i = 0; i < 256; i++)
    {
        assert(check_alphabet[i] < 2);
    }
}


// function that reads the alphabet, the score matrix and the gap
void read_alphabet(FILE *fd)
{
    int i, j;
    fscanf(fd, "%d\n", &num_symbols);

    /* check hex_mode or not (normal_mode) */
    int hex_mode;
    char c1;
    fscanf(fd, "%c", &c1);          /* dummy read */
    fscanf(fd, "%c", &c1);
    fseek(fd, 0, SEEK_SET);
    fscanf(fd, "%d\n", &hex_mode);  /* dummy read */
    hex_mode = (c1 != 0x20) ? 1 : 0;

    for (i = 1; i <= num_symbols; i++)
    {
        if (hex_mode)
        {
            int c;
            fscanf(fd, "%x ", &c);
            alphabet[i] = (unsigned char)c;

        } else {
            fscanf(fd, "%c ", &alphabet[i]);
        }
    }
    check_alphabet();

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
    long len_seq, len_seq2;
    unsigned char *seq, *seq2;
    float **matrix, res;
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

    hex_output_mode = 0;	// global
    if (argc > 3){
        if (strcmp(argv[3], "nonprintable") == 0)
	{
            hex_output_mode = 1;	// global
	}
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
        // number, name, length, sequence (numeroseq,name,longitud,sequencia)
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
        matrix = (float **)malloc(MAXLENALIGN * sizeof(float *));
        if (matrix == NULL)
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
        while ((i < MAXLENALIGN) && (estat == 0))
        {
            matrix[i] = (float *)malloc(MAXLENALIGN * sizeof(float));
            if (matrix[i] == NULL)
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
        path_matrix = (char **)malloc(MAXLENALIGN * sizeof(char *));
        if (path_matrix == NULL)
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
        while ((i < MAXLENALIGN) && (estat == 0))
        {
            path_matrix[i] = (char *)malloc(MAXLENALIGN * sizeof(char));
            if (path_matrix[i] == NULL)
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

        score_matrix = (float **)malloc(num_seqs * sizeof(float *));
        if (score_matrix == NULL)
        {
            estat = -1;
        }
        i = 0;
        estat = 0;
        while ((i < num_seqs) && (estat == 0))
        {
            score_matrix[i] = (float *)malloc(num_seqs * sizeof(float));
            if (score_matrix[i] == NULL)
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
        seq = (unsigned char *)malloc(MAXLENSEQ * sizeof(unsigned char));
        seq2 = (unsigned char *)malloc(MAXLENSEQ * sizeof(unsigned char));

        // there is enough memory for all variables
        // printf("\n 1");
        while (j < (num_seqs - 1))
        {
            len_seq = get_sequence_length(j);
            load_sequence(seq, j);
            // check_load_sequence(seq,len_seq);
            k = j + 1;
            printf("\n");
            while (k < num_seqs)
            {
                printf(".");
                len_seq2 = get_sequence_length(k);
                load_sequence(seq2, k);
                // check_load_sequence(seq2,len_seq2);
                res = similarity(matrix, seq, seq2, len_seq, len_seq2);
                //      check_path_matrix(matrix,len_seq,len_seq2);
                //      check_similarity(seq, seq2, len_seq, len_seq2);
                score_matrix[j][k] = res;
                score_matrix[k][j] = res;
                k++;
            }
            j++;
        }
        check_similarity_matrix(score_matrix);
        free(seq);
        free(seq2);
        // printf("\n 2");
        multiple_alignment(matrix, atoi(argv[2]));
        // second argument is the output format
        // printf("\n 3");
        i = 0;
        while (i < num_seqs)
        {
            free(score_matrix[i]);
            i++;
        }
        free(score_matrix);

        i = 0;
        while (i < MAXLENALIGN)
        {
            free(matrix[i]);
            i++;
        }
        free(matrix);

        i = 0;
        while (i < MAXLENALIGN)
        {
            free(path_matrix[i]);
            i++;
        }
        free(path_matrix);
        fclose(temp_file);
        remove(clusters_filename);
        remove(temp_file_name);
        remove(cluster_filename_1);
        remove(cluster_filename_2);
    }
}
