/***************************************************************************
    seqs_util.h  -  description
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

#ifndef ESSEM_SEQS_UTIL_H
#define ESSEM_SEQS_UTIL_H

extern int seqs_getInitialGaps(unsigned char * /*dna*/);
extern int seqs_getFinalGaps(unsigned char * /*dna */);
extern unsigned char *seqs_insertInitialGaps(unsigned char * /*dna*/, int /*number*/);
extern unsigned char *seqs_insertFinalGaps(unsigned char * /*dna*/, int /*number*/);

#endif
