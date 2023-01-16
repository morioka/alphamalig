/***************************************************************************
    seqs_util.c  -  description
    -------------------
    begin                : Sun Sep 9 2001
    copyright            : (C) 2001 by Roman Roset
    email                : roman.roset@menta.net

    www.fib.upc.es
    www.lsi.upc.es/~alggen
    www.cepba.upc.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *  This program is property of U.P.C (Technical University of Catalonia)  *
 *  This program is made in CIRI (cepba-ibm ,European Center for           *
 *  Parallelism of Barcelona) like project for the Group ALGGEN            *
 *  (Algorithmics and Genetics Group).                                     *
 *                                                                         *
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqs_util.h"
#include "util.h"

int seqs_getInitialGaps(char *sequence)
{
        int i = 0;

        while (sequence[i] == '-')
                i++;
        return i;
}

int seqs_getFinalGaps(char *sequence)
{
        int l = strlen(sequence);
        int i = l - 1;

        while (sequence[i] == '-')
                i--;
        return l - i - 1;
}

char *
seqs_insertInitialGaps(char *sequence, int n)
{
        int l, new_size;
        char *result;

        if (!n)
                return sequence;
        l = strlen(sequence);
        new_size = l + n + 1;
        result = malloc(new_size);
        memcpy(result + n, sequence, l + 1);
        memset(result, '-', n);
        FREE(sequence);
        return result;
}

char *
seqs_insertFinalGaps(char *sequence, int n)
{
        int l, new_size;
        char *result;

        if (!n)
                return sequence;
        l = strlen(sequence);
        new_size = l + n + 1;
        result = malloc(new_size);
        memcpy(result, sequence, l);
        memset(result + l, '-', n);
        result[l + n] = '\0';
        FREE(sequence);
        return result;
}
