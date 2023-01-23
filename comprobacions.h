

#ifndef _COMPROBACIONS_H_
#define _COMPROBACIONS_H_

void comprovar_fitxer_temporal();
void check_reading_alphabet();
void comprovar_alineament_optim(float res, char *seq1, char *seq2, long longseq1, int k);
void check_similarity(char *seq1, char *seq2, long longseq1, long longseq2);

void check_load_sequence(char *seq, long longseq);
void comprovar_matriu_similaritats(float **matriu);
void check_cluster_info(float **info);

void check_path_matrix(float **, int l1, int l2);

#endif