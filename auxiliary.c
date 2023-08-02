/*************************************************************************/
// File: auxiliary.c
// Author: Xavier Sol Acha
// Modified: Xavier Messeguer  07-01-03
// Description: Utility functions for the application
/*************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "align.h"
#include "const.h"
#include "auxiliary.h"
#include "pairs.h"

// Implementation of functions

// look at which position in the penalty matrix. there goes the symbol
int symbol_index(unsigned char c)
{
    int i = 1;
//        while ((alphabet[i] != toupper(c)) && (i < num_symbols))
    while ((alphabet[i] != c) && (i < num_symbols))
        i++;
//        if (alphabet[i] == toupper(c))
    if (alphabet[i] == c)
    {
        return i;
    }
    else
    {
        printf("Answer: The symbol %c is not in the alphabet", c);
        exit(1);
    }
}

/* belong_alphabet **********************************************************/
// It tells us if the character we pass to it as a parameter corresponds to some of the symbols of the alphabet
int belong_alphabet(unsigned char c)
{
    int i = 1;

//        c = toupper(c);
    while ((i < num_symbols) && (alphabet[i] != c))
        i++;
    return (alphabet[i] == c);
}

/* real_max *****************************************************************/
// Returns at most three real numbers.
float real_max(float x, float y, float z, unsigned char *c)
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

/* count_sequences ******************************************************/
// Counts the number of sequences in the file.
int count_sequences(void)
{
    unsigned char line[MAXLINE + 1];
    int n_seqs = 0;

    while (fgets(line, MAXLINE + 1, input_file) != NULL)
    {
        if (line[0] == '>')
        {
            n_seqs++;
        }
    }
    fseek(input_file, 0, SEEK_SET);
    return (n_seqs);
}

void show_matrix(float **m, int files, int cols)
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

int generate_prefix_files(char *fit_name)
{
    int val, i;
    time_t t;
    int modul;

    i = 0;
    val = 0;
    while ((i < 13) && (fit_name[i] != '\n') && (fit_name[i] != '\0'))
    {
        val = val + toascii(fit_name[i]);
        i++;
    }
    val = val + time(&t);
    modul = val % 100000;

    return (modul);
}

void fill_string_prefix(char *pref, int num_pref)
{
    int aux;

    aux = num_pref / 10000;
    pref[0] = to_character(aux);

    aux = num_pref / 1000;
    aux = aux % 10;
    pref[1] = to_character(aux);

    aux = num_pref / 100;
    aux = aux % 10;
    pref[2] = to_character(aux);

    aux = num_pref / 10;
    aux = aux % 10;
    pref[3] = to_character(aux);

    aux = num_pref;
    aux = aux % 10;
    pref[4] = to_character(aux);
    pref[5] = '\0';
}

char to_character(int num)
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
