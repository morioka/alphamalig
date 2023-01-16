/*************************************************************************/
// Fitxer: multiple.h
// Autor: Xavier Sol Acha
// Descripci: Aqu estan les capsaleres de les funcions per fer l'aliniament
// mltiple de seqncies
/*************************************************************************/

#ifndef _MULTIPLE_H_
#define _MULTIPLE_H_

// Temporary Files
extern FILE *fitxer_clusters;
extern FILE *fitxer_cluster_1;
extern FILE *fitxer_cluster_2;

extern char nom_fitxer_clusters[13];
extern char nom_fitxer_cluster_1[13];
extern char nom_fitxer_cluster_2[13];

extern int **seqs; // Modify to fix the bug

typedef struct
{
        int num_seqs;   // Number of cluster sequences
        long long_seqs; // Length of cluster sequences
        long puntuacio; // Score of cluster
} cluster;

extern cluster **clusters;            // Cluster vector
extern long *pos_seq;                 // Cont la posici de la seqencia al fitxer on guardem les seqncies semialiniades
extern long *pos_seq_fitxer_clusters; // Cont la posici de la seqncia al fitxer del cluster que ser utilitzat per aliniar
extern int *seq_cluster;              // Cont a quin cluster pertany cada seqncia
extern int *cluster_equivalent;       // Cont equivalncies entre clusters

float calcul_simil_ini(float *freq, int n);
float calcul_simil(float *freq1, float *freq2);
float intracluster(float *freq);
float intercluster(float *freq1, float *freq2);

// Definition of prototypes
void alineament_multiple(float **matriu, int output_format);
void trobar_clusters_propers(int *i, int *j);
void alineament_clusters(float **matriu, int i, int j, float **info1, float **info2, char *linia_orig, char *linia_desti, char *resultat_cluster1, char *resultat_cluster2, long *puntuacio);
long afegir_sequencia_cluster(FILE *fitxer, int num_seq);
void crear_fitxers_clusters(int i, int j);
void similitud_clusters(float **matriu, long long_cluster1, long long_cluster2, float **info1, float **info2, long *puntuacio,
                        int, int);
void construir_info_cluster(int num_cluster, float **info, FILE *f);
void alinea_clusters(int i, int j, char *res1, char *res2, int *len);
void reconstruir_nou_cluster(char *linia_orig, char *linia_desti, char *resultat_cluster1, char *resultat_cluster2, int i, int j,
                             int len);
void recalcular_matriu(int i, int j);

#endif
