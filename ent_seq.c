/*************************************************************************/
// File: ent_seq.c
// Author: Xavier Sol Acha
// Description: Here are the functions that will read 
// and process the DNA strings found in the input file
/*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "align.h"
#include "ent_seq.h"
#include "const.h"
#include "auxiliar.h"

// Implementation of the functions

/* read_sequences_file **********************************************/
// This function makes as many calls to load_sequence as there are 
// sequences in the original file. Returns the number of sequences in the file.
int read_sequences_file(void)
{
        int i, estat = 1;

        num_seqs = count_sequences(); // Leave the pointer at the beginning of the file
        if ((num_seqs <= MAXSEQ) && (num_seqs > 1))
        {
                i = 0;
                while ((i < num_seqs) && (estat == 1))
                {
                        estat = write_sequence_tmp(i);

                        i++;
                }
        }
        else
        {
                num_seqs = -1;
        }
        return (num_seqs);
}

/* write_sequence_tmp ******************************************************/
// It reads the following sequence that appears in the file 
// and writes it to the temporary file in my format
int write_sequence_tmp(int num_seq)
{
        static char line[MAXLONGSEQ + 1];
        static char seq_nucl[MAXLONGSEQ + 1];
        static char nom[MAXLONGNOM + 1];
        char c;
        int i, long_seq = 0, long_nom = 0, estat = 1;

        while (((*line) != '>') && !(feof(input_file)))
        {
                fgets(line, MAXLINE, input_file);
        }

        strncpy(nom, line + 1, MAXLONGNOM);
        nom[MAXLONGNOM] = EOS;

        while (fgets(line, MAXLINE + 1, input_file) && (long_seq < MAXLONGSEQ) && (*line != '>'))
        {
                i = 0;
                c = line[i];
                while ((i <= MAXLINE) && (c != '\n') && (c != EOS) && (c != '>'))
                {
//                        c = toupper(line[i]);
                        c = line[i];
                        if (belong_alphabet(c))
                        {
//                                seq_nucl[long_seq] = toupper(c);
                                seq_nucl[long_seq] = c;
                                long_seq++;
                        }
                        i++;
                }
        }
        if (long_seq == MAXLONGSEQ)
        {
                printf("WARNING: Sequence %s may have been cut up to %d nucleotides\n", nom, MAXLONGSEQ);
        }
        seq_nucl[long_seq] = EOS;

        fprintf(temp_file, ">%d\n", num_seq);
        fprintf(temp_file, "%s\n", nom);
        fprintf(temp_file, "@%d\n", long_seq);
        fprintf(temp_file, "#%s\n", seq_nucl);

        return (estat);
}

/* dona_longitud_sequencia ******************************************************/
// Returns the length of the sequence "num_seq", reading it from the temporary file
long dona_longitud_seq(int num_seq)
{
        static char line[MAXLINE + 1];
        static char aux[10];
        int i = 1, trobat = FALSE;

        fseek(temp_file, 0, SEEK_SET);
        while (trobat == FALSE)
        {
                while ((*line) != '>')
                {
                        fgets(line, MAXLINE, temp_file);
                }
                strcpy(aux, line + 1);
                if (atoi(aux) == num_seq)
                {
                        trobat = TRUE;
                }
                else
                {
                        fgets(line, MAXLINE, temp_file);
                }
        }
        while ((*line) != '@')
        {
                fgets(line, MAXLINE, temp_file);
        }
        strcpy(aux, line + 1);

        return (atol(aux));
}

/* carregar_sequencia ******************************************************/
// Returns the encoding of the sequence "num_seq" in the array "seq"
void carregar_sequencia(char *seq, int num_seq)
{
        static char line[MAXLONGSEQ + 1];
        static char aux[10];
        int i = 1, trobat = FALSE;

        fseek(temp_file, 0, SEEK_SET);
        while (trobat == FALSE)
        {
                while ((*line) != '>')
                {
                        fgets(line, MAXLONGSEQ, temp_file);
                }
                strcpy(aux, line + 1);
                if (atoi(aux) == num_seq)
                {
                        trobat = TRUE;
                }
                else
                {
                        fgets(line, MAXLONGSEQ, temp_file);
                }
        }
        while ((*line) != '#')
        {
                fgets(line, MAXLONGSEQ, temp_file);
        }
        strcpy(seq, line + 1);
}

/* load_sequence_exact_position ***************************************/
// Returns the encoding of the sequence pointed to by the read pointer. 
// This pointer must point to the line before the encoding of the sequence
void load_sequence_exact_position(char *seq)
{
        static char line[MAXLONGSEQ + 1];

        fgets(line, MAXLINE, temp_file);
        strcpy(seq, line + 1);
}

/* load_sequence_name ******************************************************/
// Return the first 10 characters of the name "num_seq" to the array "nom(name)"
void load_sequence_name(char *nom, int num_seq)
{
        static char line[MAXLINE + 1];
        static char aux[10];
        int i = 0, trobat = FALSE;

        fseek(temp_file, 0, SEEK_SET);
        while (trobat == FALSE)
        {
                while ((*line) != '>')
                {
                        fgets(line, MAXLINE, temp_file);
                }
                strcpy(aux, line + 1);
                if (atoi(aux) == num_seq)
                {
                        trobat = TRUE;
                }
                else
                {
                        fgets(line, MAXLINE, temp_file);
                }
        }

        fgets(line, MAXLINE, temp_file);
        while ((i < 10) && (line[i] != EOS))
        {
                nom[i] = line[i];
                i++;
        }
        nom[i] = EOS;
}
