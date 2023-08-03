/*************************************************************************/
// File: sort_seq.c
// Author: Xavier Sol Acha
// Description: Here are the functions to write the output of the program
// in a file.
/*************************************************************************/

#include <stdio.h>
#include "sort_seq.h"
#include "align.h"
#include "auxiliary.h"
#include "ent_seq.h"
#include "const.h"
#include "multiple.h"

#ifdef	HEX_ONLY_MODE
int hex_mode = 1;
#else
int hex_mode = 0;
#endif

void write_align_output_file(int output_format)
{
    int k = 0, i, j, l;
    unsigned char asteriscs[MAXLENSEQ];
    unsigned char seq_out[MAXLENSEQ];
//    output_format = 1;      // cut sequences
    output_format = 0;      // do not cut sequences
    if (output_format == 0) /* do not cut sequences */
    {
        printf("Number of sequences=%d  Alignment length=%ld  Alignment score=%ld\n", clusters[0]->num_seqs,
               clusters[0]->len_seqs, clusters[0]->score);
        i = 0;
        while (i < clusters[0]->len_seqs)
        {
            asteriscs[i] = '*';
            i++;
        }
        asteriscs[i] = EOS;
        while (k < clusters[0]->num_seqs)
        {
            write_sequence_cluster(seqs[0][k], asteriscs);
            k++;
        }
        i = 0;
        while (i <= 10)
        {
            seq_out[i] = ' ';
            i++;
        }
        k = 0;
        while (k < clusters[0]->len_seqs)
        {
            //if ((asteriscs[k] == 'A') || (asteriscs[k] == 'C') || (asteriscs[k] == 'G') || (asteriscs[k] == 'T'))
            //{
            //    seq_out[i] = '*';
            //}
            //else
            //{
            seq_out[i] = ' ';
            //}
            i++;
            k++;
        }
        seq_out[i] = EOS;
        printf("%s\n", seq_out);
    }
    else /* cut the sequences */
    {
        printf("\nNumber of sequences=%d  Alignment length=%ld  Alignment score=%ld\n", clusters[0]->num_seqs,
               clusters[0]->len_seqs, clusters[0]->score);
        l = 0;
        while (l < clusters[0]->len_seqs)
        {
            i = 0;
            while (i < 58)
            {
                asteriscs[i] = '*';
                i++;
            }
            asteriscs[i] = EOS;
            k = 0;
            while (k < clusters[0]->num_seqs)
            {
                write_cut_cluster_sequence(seqs[0][k], asteriscs, l);
                k++;
            }
            i = 0;
            while (i <= 10)
            {
                seq_out[i] = ' ';
                i++;
            }
            k = 0;
            while (k < 58)
            {
                //if ((asteriscs[k] == 'A') || (asteriscs[k] == 'C') || (asteriscs[k] == 'G') || (asteriscs[k] == 'T'))
                //{
                //  seq_out[i] = '*';
                //}
                //else
                //{
                seq_out[i] = ' ';
                //}
                i++;
                k++;
            }
            seq_out[i] = EOS;
            printf("%s\n\n", seq_out);
            l = l + 58;
        }
    }
}

void write_cut_cluster_sequence(int num_seq, unsigned char *asteriscs, int pos)
{
    unsigned char seq[MAXLENSEQ + 1], name[15], seq_out[MAXLENSEQ + 1];
    int i = 0, j = 0, k = 0, trobat_fi = 0;

    // First we copy the name
    fseek(clusters_file, pos_seq[num_seq], SEEK_SET);
    load_sequence_name(name, num_seq);
    i = 0;
    trobat_fi = 0;
    while (i <= 10)
    {
        if ((name[i] == EOS) || (name[i] == '\n') || (trobat_fi == 1))
        {
            seq_out[i] = ' ';
            trobat_fi = 1;
        }
        else
        {
            seq_out[i] = name[i];
        }
        i++;
    }

    j = pos;
    k = 0;
    fgets(seq, MAXLENSEQ, clusters_file);
    while ((seq[j] != EOS) && (seq[j] != '\n') && (i <= MAXLENSEQ) && ((j - pos) < 58))
    {
        seq_out[i] = seq[j];
        if (asteriscs[k] == '*')
        {
            asteriscs[k] = seq[j];
        }
        else
        {
            if (asteriscs[k] != '-')
            {
                if (asteriscs[k] != seq[j])
                {
                    asteriscs[k] = '-';
                }
            }
        }
        i++;
        j++;
        k++;
    }
    seq_out[i] = EOS;
    asteriscs[k] = EOS;
    printf("%s\n", seq_out);
}

void write_sequence_cluster(int num_seq, unsigned char *asteriscs)
{
    unsigned char seq[MAXLENSEQ + 1], name[15], seq_out[MAXLENSEQ + 1];
    int i = 0, j = 0, trobat_fi = 0;

//    int hex_mode = 1;
    int hex_mode = 0;
    hex_mode = hex_output_mode;

    // first, copy the name
    fseek(clusters_file, pos_seq[num_seq], SEEK_SET);
    load_sequence_name(name, num_seq);
    i = 0;
    trobat_fi = 0;
    while (i <= 10)
    {
        if ((name[i] == EOS) || (name[i] == '\n') || (trobat_fi == 1))
        {
            seq_out[i] = ' ';
            trobat_fi = 1;
        }
        else
        {
            seq_out[i] = name[i];
        }
        i++;
    }

    j = 0;
    fgets(seq, MAXLENSEQ, clusters_file);
    while ((seq[j] != EOS) && (seq[j] != '\n') && (i <= MAXLENSEQ))
    {
        seq_out[i] = seq[j];
        if (asteriscs[j] == '*')
        {
            asteriscs[j] = seq[j];
	    if (hex_mode) {
                unsigned char c, c0, c1;
		c = seq_out[i];
		c0 = c % 16;
		c1 = c / 16;
		if (c0 > 9)
		{
		    c0 = c0 - 10  + 'a';
	        } else {
		    c0 = c0 + '0';
		}
		if (c1 > 9)
		{
		    c1 = c1 - 10  + 'a';
	        } else {
		    c1 = c1 + '0';
		}
		if (c == 0x2d)
		{
		    c0 = '-';
		    c1 = '-';
		}
		seq_out[i] = c1;
		seq_out[i+1] = c0;
		seq_out[i+2] = ' ';
                i++; i++;
            }
        }
        else
        {
            if (asteriscs[j] != '-')
            {
                if (asteriscs[j] != seq[j])
                {
                    asteriscs[j] = '-';
                }
            }
	    if (hex_mode) {
                unsigned char c, c0, c1;
		c = seq_out[i];
		c0 = c % 16;
		c1 = c / 16;
		if (c0 > 9)
		{
		    c0 = c0 - 10  + 'a';
	        } else {
		    c0 = c0 + '0';
		}
		if (c1 > 9)
		{
		    c1 = c1 - 10  + 'a';
	        } else {
		    c1 = c1 + '0';
		}
		if (c == 0x2d)
		{
		    c0 = '-';
		    c1 = '-';
		}
		seq_out[i] = c1;
		seq_out[i+1] = c0;
		seq_out[i+2] = ' ';
                i++; i++;
            }
        }
        i++;
        j++;
    }
    seq_out[i] = EOS;
    asteriscs[i] = EOS;
    printf("%s\n", seq_out);
}
