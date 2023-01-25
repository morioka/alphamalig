

#ifndef _CHECKS_H_
#define _CHECKS_H_

void check_temp_file();
void check_reading_alphabet();
void check_optimal_alignment(float res, char *seq1, char *seq2, long lenseq1, int k);
void check_similarity(char *seq1, char *seq2, long lenseq1, long lenseq2);

void check_load_sequence(char *seq, long lenseq);
void check_similarity_matrix(float **matrix);
void check_cluster_info(float **info);

void check_path_matrix(float **, int l1, int l2);

#endif