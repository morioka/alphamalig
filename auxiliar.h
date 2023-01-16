/*************************************************************************/
// Fitxer: auxiliar.h
// Autor: Xavier Sol� Acha
// Descripci�: Fitxer de capsalera per les funcions d'utilitat per
// l'aplicaci�.
/*************************************************************************/

// Definici� dels prototips

#ifndef _AUXILIAR_H_
#define _AUXILIAR_H_

int pertany_alfabet(char c);
float maxim_real(float x, float y, float z, char *c);
int compta_sequencies(void);
void mostrar_matriu(float **m, int files, int cols); // funcio chivato ********************************
int generar_prefix_fitxers(char *nom_fit);
void omplir_string_prefix(char *pref, int num_pref);
char a_caracter(int num);
int indexsimbol(char c);

#endif
