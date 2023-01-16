/*************************************************************************/
// Fitxer: align.c
// Autor: Xavier Sol Acha
// Descripci: Aquest fitxer contindr el programa principal de l'aplicaci
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

FILE *fitxer_entrada;  // file of entry of the sequencies
FILE *fitxer_temporal; // temporary file that stores the sequencies
FILE *fitxer_alfabet;  // File that stores the alphabet in the appropriate format

char nom_fitxer_temporal[13]; // Temporary file name

int numsimb;              // number of alphabet symbols
char alfabet[MAXLONGALF]; // alphabet in positions 1,2,..,numsimb
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
void leer_alfabeto(FILE *fd)
{
    int i, j;
    fscanf(fd, "%d\n", &numsimb);

    for (i = 1; i <= numsimb; i++)
    {
        fscanf(fd, "%c ", &alfabet[i]);
    }
    //for (i = 1; i < numsimb; i++)
    //{
    //    alfabet[i] = toupper(alfabet[i]);
    //}
    for (i = 1; i <= numsimb; i++)
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
    if ((fitxer_alfabet = fopen(argv[1], "r")) == NULL)
    {
        printf("Error reading the alphabet file\n");
        exit(1);
    }
    leer_alfabeto(fitxer_alfabet);
    comprovarlecturaalfabet();
    // initialize the global variables alphabet, matpenal, numsimb
    if ((fitxer_entrada = fopen(argv[2], "r")) == NULL)
    {
        printf("Error reading the sequence file\n");
        exit(1);
    }

    // We delete temporary files. We close the ones we don't need at the moment
    num_pref = generar_prefix_fitxers(argv[1]);
    omplir_string_prefix(pref, num_pref);
    sprintf(nom_fitxer_temporal, "%s1.tmp", pref);
    sprintf(nom_fitxer_clusters, "%s2.tmp", pref);
    sprintf(nom_fitxer_cluster_1, "%s3.tmp", pref);
    sprintf(nom_fitxer_cluster_2, "%s4.tmp", pref);
    fitxer_temporal = fopen(nom_fitxer_cluster_1, "w");
    fclose(fitxer_temporal);
    fitxer_temporal = fopen(nom_fitxer_cluster_2, "w");
    fclose(fitxer_temporal);
    fitxer_temporal = fopen(nom_fitxer_clusters, "w");
    fclose(fitxer_temporal);
    fitxer_temporal = fopen(nom_fitxer_temporal, "w");
    if (llegir_sequencies_fitxer() == -1)
    // have the sequences in the temporary file in a special format
    // number, name, length, sequence (numeroseq,nom,longitud,sequencia)
    {
        printf("Too many (>%d) or too few (<2) input sequences\n", MAXSEQ);
        remove(nom_fitxer_clusters);
        remove(nom_fitxer_cluster_1);
        remove(nom_fitxer_cluster_2);
        remove(nom_fitxer_temporal);
        exit(-1);
    }
    else // correct number of sequences
    {

        fclose(fitxer_temporal);
        //       comprovar_fitxer_temporal();
        fitxer_temporal = fopen(nom_fitxer_temporal, "r");

        estat = 0;
        matriu = (float **)malloc(MAXLONGALIN * sizeof(float *));
        if (matriu == NULL)
        {
          estat = -1;
          printf("Out of memory\n");
          remove(nom_fitxer_clusters);
          remove(nom_fitxer_cluster_1);
          remove(nom_fitxer_cluster_2);
          remove(nom_fitxer_temporal);
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
            remove(nom_fitxer_clusters);
            remove(nom_fitxer_cluster_1);
            remove(nom_fitxer_cluster_2);
            remove(nom_fitxer_temporal);
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
          remove(nom_fitxer_clusters);
          remove(nom_fitxer_cluster_1);
          remove(nom_fitxer_cluster_2);
          remove(nom_fitxer_temporal);
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
            remove(nom_fitxer_clusters);
            remove(nom_fitxer_cluster_1);
            remove(nom_fitxer_cluster_2);
            remove(nom_fitxer_temporal);
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
              remove(nom_fitxer_clusters);
              remove(nom_fitxer_cluster_1);
              remove(nom_fitxer_cluster_2);
              remove(nom_fitxer_temporal);
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
            // comprovar_carregar_sequencia(seq,long_seq);
            k = j + 1;
            printf("\n");
            while (k < num_seqs)
            {
                printf(".");
                long_seq2 = dona_longitud_seq(k);
                carregar_sequencia(seq2, k);
                // comprovar_carregar_sequencia(seq2,long_seq2);
                res = similitud(matriu, seq, seq2, long_seq, long_seq2);
                //      comprova_matriu_cami(matriu,long_seq,long_seq2);
                //      comprovar_similitud(seq, seq2, long_seq, long_seq2);
                matriu_puntuacions[j][k] = res;
                matriu_puntuacions[k][j] = res;
                k++;
            }
            j++;
        }
        comprovar_matriu_similaritats(matriu_puntuacions);
        free(seq);
        free(seq2);
        // printf("\n 2");
        alineament_multiple(matriu, atoi(argv[2]));
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
        fclose(fitxer_temporal);
        remove(nom_fitxer_clusters);
        remove(nom_fitxer_temporal);
        remove(nom_fitxer_cluster_1);
        remove(nom_fitxer_cluster_2);
    }
}
