/*************************************************************************/
// File: sort_seq.c
// Author: Xavier Sol Acha
// Description: Here are the headers of the functions to write the
// program output to a file.
/*************************************************************************/

#include "parells.h"

#ifndef _SORT_SEQ_H_
#define _SORT_SEQ_H_

void write_align_output_file(int output_format);
void write_sequence_cluster(int num_seq, char *asteriscs);
void write_cut_cluster_sequence(int num_seq, char *asteriscs, int pos);

#endif
