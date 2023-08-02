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
#include "auxiliary.h"

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
    static unsigned char line[MAXLENSEQ + 1];
    static unsigned char seq_nucl[MAXLENSEQ + 1];
    static char name[MAXLENNAME + 1];
    unsigned char c;
    int i, len_seq = 0, len_name = 0, estat = 1;
    unsigned char c0, c1;

    while (((*line) != '>') && !(feof(input_file)))
    {
        fgets(line, MAXLINE, input_file);
    }

    strncpy(name, line + 1, MAXLENNAME);
    name[MAXLENNAME] = EOS;

    while (fgets(line, MAXLINE + 1, input_file) && (len_seq < MAXLENSEQ) && (*line != '>'))
    {
	int hex_mode = 0;
        i = 0;
        c = line[i];
	c0= line[i+1];	// dummy
	c1= line[i+2];	// if hex_mode, this value is 0x20.
	hex_mode = (c1 == 0x20) ? 1 : 0;

	if (hex_mode)
	{
	    printf("hex_mode\n");
	}

        while ((i <= MAXLINE) && (c != '\n') && (c != EOS) && (c != '>'))
        {
            //c = toupper(line[i]);
            c = line[i];
	    if (hex_mode)
            {
                c0 = line[i+1];
		c1 = line[i+2]; // dummy
		
		if (c > 0x61)		// a,b,c,d,e,f
		{
		    c -= 0x20;
		}

		if (c > 0x41)
	        {
	            c = c - 0x41 + 10;	// A,B,C,D,E,F
		}
	        else
		{
		    c = c - 0x30;	// 0-9
		}
					//
		if (c0> 0x61)		// a,b,c,d,e,f
		{
	            c0 -= 0x20;
		}

		if (c0> 0x41)
	        {
	            c0= c0 - 0x41 + 10;	// A,B,C,D,E,F
		}
	        else
		{
		    c0= c0 - 0x30;	// 0-9
		}

		c = (char)((c << 4) + c0);
	    }

            if (belong_alphabet(c))
            {
                //seq_nucl[len_seq] = toupper(c);
                seq_nucl[len_seq] = c;
                len_seq++;
            }
            i++;
	    if (hex_mode)
            {
                c = c1;
                i++;
                i++;
	    }
	     
        }
    }
    if (len_seq == MAXLENSEQ)
    {
        printf("WARNING: Sequence %s may have been cut up to %d nucleotides\n", name, MAXLENSEQ);
    }
    seq_nucl[len_seq] = EOS;

    fprintf(temp_file, ">%d\n", num_seq);
    fprintf(temp_file, "%s\n", name);
    fprintf(temp_file, "@%d\n", len_seq);
    fprintf(temp_file, "#%s\n", seq_nucl);

    return (estat);
}

/* get_sequence_length ******************************************************/
// Returns the length of the sequence "num_seq", reading it from the temporary file
long get_sequence_length(int num_seq)
{
    static unsigned char line[MAXLINE + 1];
    static unsigned char aux[10];
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

/* load_sequence ******************************************************/
// Returns the encoding of the sequence "num_seq" in the array "seq"
void load_sequence(unsigned char *seq, int num_seq)
{
    static unsigned char line[MAXLENSEQ + 1];
    static unsigned char aux[10];
    int i = 1, trobat = FALSE;

    fseek(temp_file, 0, SEEK_SET);
    while (trobat == FALSE)
    {
        while ((*line) != '>')
        {
            fgets(line, MAXLENSEQ, temp_file);
        }
        strcpy(aux, line + 1);
        if (atoi(aux) == num_seq)
        {
            trobat = TRUE;
        }
        else
        {
            fgets(line, MAXLENSEQ, temp_file);
        }
    }
    while ((*line) != '#')
    {
        fgets(line, MAXLENSEQ, temp_file);
    }
    strcpy(seq, line + 1);
}

/* load_sequence_exact_position ***************************************/
// Returns the encoding of the sequence pointed to by the read pointer.
// This pointer must point to the line before the encoding of the sequence
void load_sequence_exact_position(unsigned char *seq)
{
    static unsigned char line[MAXLENSEQ + 1];

    fgets(line, MAXLINE, temp_file);
    strcpy(seq, line + 1);
}

/* load_sequence_name ******************************************************/
// Return the first 10 characters of the name "num_seq" to the array "name"
void load_sequence_name(unsigned char *name, int num_seq)
{
    static unsigned char line[MAXLINE + 1];
    static unsigned char aux[10];
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
        name[i] = line[i];
        i++;
    }
    name[i] = EOS;
}
