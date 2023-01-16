/*************************************************************************/
// Fitxer: ent_seq.c
// Autor: Xavier Sol Acha
// Descripci: Aqu estan les funcions que llegiran i processaran les
// cadenes d'ADN que es trobin al fitxer d'entrada
/*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "align.h"
#include "ent_seq.h"
#include "const.h"
#include "auxiliar.h"

// Implementaci de les funcions

/* llegir_sequencies_fitxer **********************************************/
// Aquesta funci fa tantes crides a carregar_sequencia com seqncies hi
// hagi en el fitxer original.Retorna el nombre de seqncies que hi ha al
// fitxer.
int llegir_sequencies_fitxer(void)
{
        int i, estat = 1;

        num_seqs = compta_sequencies(); // Deixa el punter a l'inici del fitxer
        if ((num_seqs <= MAXSEQ) && (num_seqs > 1))
        {
                i = 0;
                while ((i < num_seqs) && (estat == 1))
                {
                        estat = escriure_sequencia_tmp(i);

                        i++;
                }
        }
        else
        {
                num_seqs = -1;
        }
        return (num_seqs);
}

/* escriure_sequencia_tmp ******************************************************/
// Llegeix la segent seqncia que apareix al fitxer i l'escriu al fitxer temporal
// en el meu format
int escriure_sequencia_tmp(int num_seq)
{
        static char linia[MAXLONGSEQ + 1];
        static char seq_nucl[MAXLONGSEQ + 1];
        static char nom[MAXLONGNOM + 1];
        char c;
        int i, long_seq = 0, long_nom = 0, estat = 1;

        while (((*linia) != '>') && !(feof(fitxer_entrada)))
        {
                fgets(linia, MAXLINIA, fitxer_entrada);
        }

        strncpy(nom, linia + 1, MAXLONGNOM);
        nom[MAXLONGNOM] = EOS;

        while (fgets(linia, MAXLINIA + 1, fitxer_entrada) && (long_seq < MAXLONGSEQ) && (*linia != '>'))
        {
                i = 0;
                c = linia[i];
                while ((i <= MAXLINIA) && (c != '\n') && (c != EOS) && (c != '>'))
                {
                        c = toupper(linia[i]);
                        if (pertany_alfabet(c))
                        {
                                seq_nucl[long_seq] = toupper(c);
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

        fprintf(fitxer_temporal, ">%d\n", num_seq);
        fprintf(fitxer_temporal, "%s\n", nom);
        fprintf(fitxer_temporal, "@%d\n", long_seq);
        fprintf(fitxer_temporal, "#%s\n", seq_nucl);

        return (estat);
}

/* dona_longitud_sequencia ******************************************************/
// Retorna la longitud de la sequencia "num_seq", llegint-la del fitxer temporal
long dona_longitud_seq(int num_seq)
{
        static char linia[MAXLINIA + 1];
        static char aux[10];
        int i = 1, trobat = FALSE;

        fseek(fitxer_temporal, 0, SEEK_SET);
        while (trobat == FALSE)
        {
                while ((*linia) != '>')
                {
                        fgets(linia, MAXLINIA, fitxer_temporal);
                }
                strcpy(aux, linia + 1);
                if (atoi(aux) == num_seq)
                {
                        trobat = TRUE;
                }
                else
                {
                        fgets(linia, MAXLINIA, fitxer_temporal);
                }
        }
        while ((*linia) != '@')
        {
                fgets(linia, MAXLINIA, fitxer_temporal);
        }
        strcpy(aux, linia + 1);

        return (atol(aux));
}

/* carregar_sequencia ******************************************************/
// Retorna la codificaci de la sequencia "num_seq" a l'array "seq"
void carregar_sequencia(char *seq, int num_seq)
{
        static char linia[MAXLONGSEQ + 1];
        static char aux[10];
        int i = 1, trobat = FALSE;

        fseek(fitxer_temporal, 0, SEEK_SET);
        while (trobat == FALSE)
        {
                while ((*linia) != '>')
                {
                        fgets(linia, MAXLONGSEQ, fitxer_temporal);
                }
                strcpy(aux, linia + 1);
                if (atoi(aux) == num_seq)
                {
                        trobat = TRUE;
                }
                else
                {
                        fgets(linia, MAXLONGSEQ, fitxer_temporal);
                }
        }
        while ((*linia) != '#')
        {
                fgets(linia, MAXLONGSEQ, fitxer_temporal);
        }
        strcpy(seq, linia + 1);
}

/* carregar_sequencia_posicio_exacte ***************************************/
// Retorna la codificaci de la sequencia on apunta el punter de lectura. Aquest
// punter ha d'apuntar a la lnia anterior a la codificaci de la seqncia
void carregar_sequencia_posicio_exacte(char *seq)
{
        static char linia[MAXLONGSEQ + 1];

        fgets(linia, MAXLINIA, fitxer_temporal);
        strcpy(seq, linia + 1);
}

/* carregar_nom_sequencia ******************************************************/
// Return the first 10 characters of the name "num_seq" to the array "nom(name)"
void carregar_nom_sequencia(char *nom, int num_seq)
{
        static char linia[MAXLINIA + 1];
        static char aux[10];
        int i = 0, trobat = FALSE;

        fseek(fitxer_temporal, 0, SEEK_SET);
        while (trobat == FALSE)
        {
                while ((*linia) != '>')
                {
                        fgets(linia, MAXLINIA, fitxer_temporal);
                }
                strcpy(aux, linia + 1);
                if (atoi(aux) == num_seq)
                {
                        trobat = TRUE;
                }
                else
                {
                        fgets(linia, MAXLINIA, fitxer_temporal);
                }
        }

        fgets(linia, MAXLINIA, fitxer_temporal);
        while ((i < 10) && (linia[i] != EOS))
        {
                nom[i] = linia[i];
                i++;
        }
        nom[i] = EOS;
}
