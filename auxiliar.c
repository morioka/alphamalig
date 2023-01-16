/*************************************************************************/
// Fitxer: auxiliar.c
// Autor: Xavier Sol Acha
// Modificat: Xavier Messeguer  07-01-03
// Descripci: Funcions d'utilitat per l'aplicaci
/*************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "align.h"
#include "const.h"
#include "auxiliar.h"
#include "parells.h"

// Implementaci de les funcions

// mira a quina posicio de la matriu de penal. hi va el simbol
int indexsimbol(char c)
{
        int i = 1;
        //  printf("Miro simbol %c\n",c);
        while ((alfabet[i] != toupper(c)) && (i < numsimb))
                i++;
        if (alfabet[i] == toupper(c))
        {
                /*    printf("\nResposta: index %d\n",i);*/
                return i;
        }
        else
        {
                printf("Resposta: el simbol %c no es a l'alfabet", c);
                exit(1);
        }
}

/* pertany_alfabet **********************************************************/
// Ens diu si el carcter que li passem com a parmetre correspon a alguns
// dels simbols de l'alfabet
int pertany_alfabet(char c)
{
        int i = 1;

        c = toupper(c);
        while ((i < numsimb) && (alfabet[i] != c))
                i++;
        return (alfabet[i] == c);
}

/* maxim_real *****************************************************************/
// Retorna el mxim de tres nombres reals.
float maxim_real(float x, float y, float z, char *c)
{
        float num_max;

        if (x < y)
        {
                num_max = y;
                *c = 'd';
        }
        else
        {
                num_max = x;
                *c = 'a';
        }
        if (num_max < z)
        {
                num_max = z;
                *c = 'e';
        }
        return (num_max);
}

/* compta_sequencies ******************************************************/
// Compta el nombre de seqncies que hi ha al fitxer.
int compta_sequencies(void)
{
        char linia[MAXLINIA + 1];
        int n_seqs = 0;

        while (fgets(linia, MAXLINIA + 1, fitxer_entrada) != NULL)
        {
                if (linia[0] == '>')
                {
                        n_seqs++;
                }
        }
        fseek(fitxer_entrada, 0, SEEK_SET);
        return (n_seqs);
}

void mostrar_matriu(float **m, int files, int cols)
{
        int i, j;

        for (i = 0; i < files; i++)
        {
                for (j = 0; j < cols; j++)
                {
                        printf("%f ", m[i][j]);
                }
                printf("\n");
        }
}

int generar_prefix_fitxers(char *nom_fit)
{
        int val, i;
        time_t t;
        int modul;

        i = 0;
        val = 0;
        while ((i < 13) && (nom_fit[i] != '\n') && (nom_fit[i] != '\0'))
        {
                val = val + toascii(nom_fit[i]);
                i++;
        }
        val = val + time(&t);
        modul = val % 100000;

        return (modul);
}

void omplir_string_prefix(char *pref, int num_pref)
{
        int aux;

        aux = num_pref / 10000;
        pref[0] = a_caracter(aux);

        aux = num_pref / 1000;
        aux = aux % 10;
        pref[1] = a_caracter(aux);

        aux = num_pref / 100;
        aux = aux % 10;
        pref[2] = a_caracter(aux);

        aux = num_pref / 10;
        aux = aux % 10;
        pref[3] = a_caracter(aux);

        aux = num_pref;
        aux = aux % 10;
        pref[4] = a_caracter(aux);
        pref[5] = '\0';
}

char a_caracter(int num)
{
        char val;

        switch (num)
        {
        case 0:
                val = '0';
                break;
        case 1:
                val = '1';
                break;
        case 2:
                val = '2';
                break;
        case 3:
                val = '3';
                break;
        case 4:
                val = '4';
                break;
        case 5:
                val = '5';
                break;
        case 6:
                val = '6';
                break;
        case 7:
                val = '7';
                break;
        case 8:
                val = '8';
                break;
        case 9:
                val = '9';
                break;
        }
        return (val);
}
