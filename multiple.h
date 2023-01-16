/*************************************************************************/
                                  //Fitxer: multiple.h
                                  //Autor: Xavier Solé Acha
                                  //Descripció: Aquí estan les capsaleres de les funcions per fer l'aliniament
                                  //múltiple de seqüències
                                  /*************************************************************************/

                                  //Fitxers temporals
                                  FILE *fitxer_clusters;
                                  FILE *fitxer_cluster_1;
                                  FILE *fitxer_cluster_2;

                                  char nom_fitxer_clusters[13];
                                  char nom_fitxer_cluster_1[13];
                                  char nom_fitxer_cluster_2[13];

                                  int **seqs; //Modificació per arreglar el bug

                                  typedef struct
                                  {
                                          int num_seqs;    //Nombre de seqüències del cluster
                                          long long_seqs;  //Longitud de les seqüències del cluster
                                          long puntuacio;  //Puntuacio del cluster
                                  }cluster;

                                  cluster **clusters;  //Vector de clusters
                                  long *pos_seq;  //Conté la posició de la seqüencia al fitxer on guardem les seqüències semialiniades
                                  long *pos_seq_fitxer_clusters; //Conté la posició de la seqüència al fitxer del cluster que serà utilitzat per aliniar
                                  int *seq_cluster;  //Conté a quin cluster pertany cada seqüència
                                  int *cluster_equivalent;  //Conté equivalències entre clusters

                                  float calcul_simil_ini(float * freq,int n);
                                  float calcul_simil(float *freq1, float *freq2 );
                                  float intracluster (float * freq);
                                  float intercluster (float *freq1,float *freq2);

                                  //Definicio dels prototips
                                  void alineament_multiple(float **matriu, int output_format);
                                  void trobar_clusters_propers(int *i, int *j);
                                  void alineament_clusters(float **matriu, int i, int j, float **info1, float **info2, char *linia_orig, char *linia_desti, char
                                  *resultat_cluster1, char *resultat_cluster2, long *puntuacio);
                                  long afegir_sequencia_cluster(FILE *fitxer, int num_seq);
                                  void crear_fitxers_clusters(int i, int j);
                                  void similitud_clusters(float **matriu, long long_cluster1, long long_cluster2, float **info1, float **info2, long *puntuacio,
                                  int, int);
                                  void construir_info_cluster(int num_cluster, float **info, FILE *f);
                                  void alinea_clusters(int i, int j, char *res1, char *res2, int *len);
                                  void reconstruir_nou_cluster(char *linia_orig, char *linia_desti, char *resultat_cluster1, char *resultat_cluster2, int i, int j,
                                  int len);
                                  void recalcular_matriu(int i, int j);
