/*************************************************************************/
// File: multiple.h
// Author: Xavier Sol Acha
// Description: Here are the headers of the functions to do the alignment
// multiple sequences
/*************************************************************************/

#ifndef _MULTIPLE_H_
#define _MULTIPLE_H_

// Temporary Files
extern FILE *clusters_file;
extern FILE *cluster_file_1;
extern FILE *cluster_file_2;

extern char clusters_filename[13];
extern char cluster_filename_1[13];
extern char cluster_filename_2[13];

extern int **seqs; // Modify to fix the bug

typedef struct
{
        int num_seqs;   // Number of cluster sequences
        long long_seqs; // Length of cluster sequences
        long score; // Score of cluster
} cluster;

extern cluster **clusters;            // Cluster vector
extern long *pos_seq;                 // position of the sequence in the file where we save the semi-aligned sequences
extern long *pos_seq_fitxer_clusters; // position of the sequence in the cluster file to be used for alignment
extern int *seq_cluster;              // which cluster each sequence belongs
extern int *cluster_equivalent;       // equivalences between clusters

float calc_simil_ini(float *freq, int n);
float calc_simil(float *freq1, float *freq2);
float intracluster(float *freq);
float intercluster(float *freq1, float *freq2);

// Definition of prototypes
void multiple_alignment(float **matrix, int output_format);
void find_nearest_clusters(int *i, int *j);
void alignment_clusters(float **matrix, int i, int j, float **info1, float **info2, char *line_orig, char *line_dest, char *result_cluster1, char *result_cluster2, long *score);
long add_sequence_cluster(FILE *fitxer, int num_seq);
void create_cluster_files(int i, int j);
void similarity_clusters(float **matrix, long long_cluster1, long long_cluster2, float **info1, float **info2, long *score,
                        int, int);
void build_info_cluster(int num_cluster, float **info, FILE *f);
void align_clusters(int i, int j, char *res1, char *res2, int *len);
void rebuild_new_cluster(char *line_orig, char *line_dest, char *result_cluster1, char *result_cluster2, int i, int j,
                             int len);
void recalculate_matrix(int i, int j);

#endif
