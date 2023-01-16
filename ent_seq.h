/*************************************************************************/
// Fitxer: ent_seq.h
// Autor: Xavier Sol� Acha
// Descripci�: Fitxer de capsalera per les funcions que llegiran i
// processaran les cadenes d'ADN que es trobin al fitxer d'entrada o al
// temporal
/*************************************************************************/

// Definici� dels prototips

#ifndef _ENT_SEQ_H_
#define _ENT_SEQ_H_

int llegir_sequencies_fitxer(void);
int escriure_sequencia_tmp(int num_seq);
long dona_longitud_seq(int num_seq);
void carregar_sequencia(char *seq, int num_seq);
void carregar_sequencia_posicio_exacte(char *seq);
void carregar_nom_sequencia(char *nom, int num_seq);

#endif
