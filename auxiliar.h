/*************************************************************************/
// File: auxiliar.h
// Author: Xavier Sol Acha
// Description: Header file for the utility functions for
// the application.
/*************************************************************************/

// Definici dels prototips

#ifndef _AUXILIAR_H_
#define _AUXILIAR_H_

int belong_alphabet(char c);
float maxim_real(float x, float y, float z, char *c);
int count_sequences(void);
void mostrar_matriu(float **m, int files, int cols); // Chivato function ********************************
int generate_prefix_files(char *nom_fit);
void fill_string_prefix(char *pref, int num_pref);
char to_character(int num);
int symbol_index(char c);

#endif
