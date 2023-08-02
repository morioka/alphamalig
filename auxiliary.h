/*************************************************************************/
// File: auxiliary.h
// Author: Xavier Sol Acha
// Description: Header file for the utility functions for
// the application.
/*************************************************************************/

// Definici dels prototips

#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

int belong_alphabet(unsigned char c);
float real_max(float x, float y, float z, unsigned char *c);
int count_sequences(void);
void show_matrix(float **m, int files, int cols); // Chivato function ********************************
int generate_prefix_files(char *fit_name);
void fill_string_prefix(char *pref, int num_pref);
char to_character(int num);
int symbol_index(unsigned char c);

#endif
