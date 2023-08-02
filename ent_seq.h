/*************************************************************************/
// File: ent_seq.h
// Author: Xavier Sol Acha
// Description: Header file for the functions that will read i
// will process the DNA strings found in the input file or the
// storm
/*************************************************************************/

// Definici dels prototips

#ifndef _ENT_SEQ_H_
#define _ENT_SEQ_H_

int read_sequences_file(void);
int write_sequence_tmp(int num_seq);
long get_sequence_length(int num_seq);
void load_sequence(unsigned char *seq, int num_seq);
void load_sequence_exact_position(unsigned char *seq);
void load_sequence_name(unsigned char *name, int num_seq);

#endif
