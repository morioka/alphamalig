/*************************************************************************/
// Fitxer: sort_seq.c
// Autor: Xavier Sol Acha
// Descripci: Aqu estan les funcions per escriure la sortida del programa
// en un fitxer.
/*************************************************************************/

#include <stdio.h>
#include "sort_seq.h"
#include "align.h"
#include "auxiliar.h"
#include "ent_seq.h"
#include "const.h"
#include "multiple.h"

void escriure_alineament_fitxer_sortida(int output_format)
{
    int k = 0, i, j, l;
    char asteriscs[MAXLONGSEQ];
    char seq_out[MAXLONGSEQ];
    output_format = 1;      // tallem sequencies
    if (output_format == 0) /* No tallem les sequencies */
    {
        printf("Number of sequnces=%d  Alignment length=%d  Alignment score=%d\n", clusters[0]->num_seqs,
              clusters[0]->long_seqs, clusters[0]->puntuacio);
        i = 0;
        while (i < clusters[0]->long_seqs)
        {
            asteriscs[i] = '*';
            i++;
        }
        asteriscs[i] = EOS;
        while (k < clusters[0]->num_seqs)
        {
            escriure_sequencia_cluster(seqs[0][k], asteriscs);
            k++;
        }
        i = 0;
        while (i <= 10)
        {
            seq_out[i] = ' ';
            i++;
        }
        k = 0;
        while (k < clusters[0]->long_seqs)
        {
            if ((asteriscs[k] == 'A') || (asteriscs[k] == 'C') || (asteriscs[k] == 'G') || (asteriscs[k] == 'T'))
            {
                seq_out[i] = '*';
            }
            else
            {
                seq_out[i] = ' ';
            }
            i++;
            k++;
        }
        seq_out[i] = EOS;
        printf("%s\n", seq_out);
    }
    else /* Tallem les sequencies */
    {
        printf("\nNumber of sequnces=%d  Alignment length=%d  Alignment score=%d\n", clusters[0]->num_seqs,
              clusters[0]->long_seqs, clusters[0]->puntuacio);
        l = 0;
        while (l < clusters[0]->long_seqs)
        {
            i = 0;
            while (i < 58)
            {
                asteriscs[i] = '*';
                i++;
            }
            asteriscs[i] = EOS;
            k = 0;
            // printf("\n 2.2.1\n");
            while (k < clusters[0]->num_seqs)
            {
                escriure_sequencia_cluster_tallada(seqs[0][k], asteriscs, l);
                k++;
            }
            i = 0;
            while (i <= 10)
            {
                seq_out[i] = ' ';
                i++;
            }
            k = 0;
            // printf("\n 2.2.2");
            while (k < 58)
            {
                if ((asteriscs[k] == 'A') || (asteriscs[k] == 'C') || (asteriscs[k] == 'G') || (asteriscs[k] == 'T'))
                {
                  seq_out[i] = '*';
                }
                else
                {
                  seq_out[i] = ' ';
                }
                i++;
                k++;
            }
            // printf("\n 2.2.3");
            seq_out[i] = EOS;
            printf("%s\n\n", seq_out);
            l = l + 58;
        }
    }

    printf("\n 2.2.4");
}

void escriure_sequencia_cluster_tallada(int num_seq, char *asteriscs, int pos)
{
    char seq[MAXLONGSEQ + 1], nom[15], seq_out[MAXLONGSEQ + 1];
    int i = 0, j = 0, k = 0, trobat_fi = 0;

    // Primer copiem el nom
    fseek(fitxer_clusters, pos_seq[num_seq], SEEK_SET);
    carregar_nom_sequencia(nom, num_seq);
    i = 0;
    trobat_fi = 0;
    while (i <= 10)
    {
        if ((nom[i] == EOS) || (nom[i] == '\n') || (trobat_fi == 1))
        {
            seq_out[i] = ' ';
            trobat_fi = 1;
        }
        else
        {
            seq_out[i] = nom[i];
        }
        i++;
    }

    j = pos;
    k = 0;
    fgets(seq, MAXLONGSEQ, fitxer_clusters);
    while ((seq[j] != EOS) && (seq[j] != '\n') && (i <= MAXLONGSEQ) && ((j - pos) < 58))
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

void escriure_sequencia_cluster(int num_seq, char *asteriscs)
{
    char seq[MAXLONGSEQ + 1], nom[15], seq_out[MAXLONGSEQ + 1];
    int i = 0, j = 0, trobat_fi = 0;

    // Primer copiem el nom
    fseek(fitxer_clusters, pos_seq[num_seq], SEEK_SET);
    carregar_nom_sequencia(nom, num_seq);
    i = 0;
    trobat_fi = 0;
    while (i <= 10)
    {
        if ((nom[i] == EOS) || (nom[i] == '\n') || (trobat_fi == 1))
        {
            seq_out[i] = ' ';
            trobat_fi = 1;
        }
        else
        {
            seq_out[i] = nom[i];
        }
        i++;
    }

    j = 0;
    fgets(seq, MAXLONGSEQ, fitxer_clusters);
    while ((seq[j] != EOS) && (seq[j] != '\n') && (i <= MAXLONGSEQ))
    {
        seq_out[i] = seq[j];
        if (asteriscs[j] == '*')
        {
            asteriscs[j] = seq[j];
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
        }
        i++;
        j++;
    }
    seq_out[i] = EOS;
    asteriscs[i] = EOS;
    printf("%s\n", seq_out);
}
