/************************************************************************/
// File: multiple.c
// Author: Xavier Sol Acha
// Description: Here is the implementation of the functions to do the alignment
// multiple sequences
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "align.h"
#include "ent_seq.h"
#include "const.h"
#include "auxiliar.h"
#include "parells.h"
#include "malloc.h"
#include "ent_seq.h"
#include "multiple.h"
#include "sort_seq.h"

// Fitxers temporals
FILE *clusters_file;
FILE *cluster_file_1;
FILE *cluster_file_2;

char clusters_filename[13];
char cluster_filename_1[13];
char cluster_filename_2[13];

int **seqs; // Modificaci per arreglar el bug

cluster **clusters;            // Vector de clusters
long *pos_seq;                 // Cont la posici de la seqencia al fitxer on guardem les seqncies semialiniades
long *pos_seq_clusters_file; // Cont la posici de la seqncia al fitxer del cluster que ser utilitzat per aliniar
int *seq_cluster;              // Cont a quin cluster pertany cada seqncia
int *cluster_equivalent;       // Cont equivalncies entre clusters

/* multiple_alignment */
/*Funci que portar a terme laliniament mltiple de les seqncies a
  //partir de la matriu de puntuacions dos a dos*/
void multiple_alignment(float **matriu, int output_format)
{
  int i = 0, j = 0, num_clusters = 0, l = 0, mem_ok = 1;
  char *line_orig, *line_dest;
  char *resultat_cluster1, *resultat_cluster2;
  float **info1, **info2;

  // Inicialitzacions
  line_orig = (char *)malloc((MAXLONGSEQ + 1) * sizeof(char));
  line_dest = (char *)malloc((MAXLONGSEQ + 1) * sizeof(char));

  // Vector que cont quines seqncies pertanyen a cada cluster
  seqs = (int **)malloc(num_seqs * sizeof(int *));

  // Vector de clusters
  clusters = (cluster **)malloc(num_seqs * sizeof(cluster *));
  if (clusters == NULL)
    mem_ok = 0;

  for (i = 0; (i < num_seqs) && (mem_ok == 1); i++)
  {
    clusters[i] = (cluster *)malloc(sizeof(cluster));
    if (clusters[i] == NULL)
    {
      mem_ok = 0;
    }
    else
    {
      clusters[i]->puntuacio = 0;
      clusters[i]->num_seqs = 1;
      clusters[i]->long_seqs = 0; // De moment no ho posem
      seqs[i] = (int *)malloc(1 * sizeof(int));
      if (seqs[i] == NULL)
      {
        mem_ok = 0;
      }
      else
      {
        seqs[i][0] = i;
      }
    }
  }

  // Vector que cont la posicio de les seqencies semi-aliniades
  if (mem_ok == 1)
  {
    pos_seq = (long *)malloc(num_seqs * sizeof(long));
    if (pos_seq == NULL)
    {
      mem_ok = 0;
    }
    else
    {
      for (i = 0; i < num_seqs; i++)
        pos_seq[i] = -1;
    }
  }

  // Vector que cont la posic de les seqncies al cluster_file_1 o cluster_file_2
  if (mem_ok == 1)
  {
    pos_seq_clusters_file = (long *)calloc(num_seqs, sizeof(long));
    if (pos_seq_clusters_file == NULL)
    {
      mem_ok = 0;
    }
  }

  // Vector que cont a quin cluster pertany cada seqencia
  if (mem_ok == 1)
  {
    seq_cluster = (int *)malloc(num_seqs * sizeof(int));
    if (seq_cluster == NULL)
    {
      mem_ok = 0;
    }
    else
    {
      for (i = 0; i < num_seqs; i++)
        seq_cluster[i] = i;
    }
  }

  // Vector que cont les equivalncies entre clusters
  if (mem_ok == 1)
  {
    // el cluster te el numero de la sequencia mes petita que conte
    cluster_equivalent = (int *)malloc(num_seqs * sizeof(int));
    if (cluster_equivalent == NULL)
    {
      mem_ok = 0;
    }
    else
    {
      for (i = 0; i < num_seqs; i++)
        cluster_equivalent[i] = i;
    }
  }

  resultat_cluster1 = (char *)malloc(MAXLONGSEQ * sizeof(char));
  if (resultat_cluster1 == NULL)
  {
    mem_ok = 0;
  }
  resultat_cluster2 = (char *)malloc(MAXLONGSEQ * sizeof(char));
  if (resultat_cluster2 == NULL)
  {
    mem_ok = 0;
  }

  if (mem_ok == 1)
  {
    info1 = (float **)malloc(MAXLONGSEQ * sizeof(float));
    if (info1 == NULL)
    {
      mem_ok = 0;
    }
    else
    {
      for (i = 0; i < MAXLONGSEQ; i++)
      {
        info1[i] = (float *)malloc(num_symbols * sizeof(float));
        if (info1[i] == NULL)
          mem_ok = 0;
      }
    }
  }

  if (mem_ok == 1)
  {
    info2 = (float **)malloc(MAXLONGSEQ * sizeof(float));
    if (info2 == NULL)
    {
      mem_ok = 0;
    }
    else
    {
      for (i = 0; i < MAXLONGSEQ; i++)
      {
        info2[i] = (float *)malloc(num_symbols * sizeof(float));
        if (info2[i] == NULL)
          mem_ok = 0;
      }
    }
  }

  num_clusters = num_seqs;

  // Comencem lalineament
  // printf("\n 2.1\n");
  if (mem_ok == 1)
  {
    while (num_clusters > 1)
    {
      printf(".");

      find_nearest_clusters(&i, &j);

      if ((clusters[i]->long_seqs >= MAXLONGSEQ) || (clusters[j]->long_seqs >= MAXLONGSEQ))
      {
        printf("Too long resulting alignment (>=2000)\n");
        exit(-1);
      }
      create_cluster_files(i, j);
      cluster_file_1 = fopen(cluster_filename_1, "r");
      cluster_file_2 = fopen(cluster_filename_2, "r");

      build_info_cluster(i, info1, cluster_file_1);
      //      check_cluster_info(info1);///////////////////
      build_info_cluster(j, info2, cluster_file_2);
      //      check_cluster_info(info2);///////////////////
      alignment_clusters(matriu, i, j, info1, info2, line_orig, line_dest, resultat_cluster1, resultat_cluster2, &clusters[i]->puntuacio);
      //     check_path_matrix(matriu,clusters[i]->long_seqs,clusters[j]->long_seqs);
      recalculate_matrix(i, j); // tambe s'ha fet dins aliniament_clusters
      fclose(cluster_file_1);
      fclose(cluster_file_2);
      num_clusters--;
    }
    clusters_file = fopen(clusters_filename, "r");
    write_align_output_file(output_format);
    fclose(clusters_file);
  }
  else
  {
    printf("Out of memory\n");
    remove(clusters_filename);
    remove(cluster_filename_1);
    remove(cluster_filename_2);
    remove(temp_file_name);
    exit(-1);
  }
}

/*build_info_cluster */
// Given a file with a cluster, populate the number of elements
// of the alphabet in each column of the cluster
void build_info_cluster(int num_cluster, float **info, FILE *f)
{
  char line[MAXLONGSEQ + 1];
  int i, j, k;

  // Initialize the structures
  j = 0;

  while (j < clusters[num_cluster]->long_seqs)
  {
    for (k = 0; k < num_symbols; k++)
      info[j][k] = 0;
    j++;
  }
  i = 0;
  fseek(f, 0, SEEK_SET);
  while (i < clusters[num_cluster]->num_seqs)
  {
    fgets(line, MAXLONGSEQ, f);
    j = 0;
    while (j < clusters[num_cluster]->long_seqs)
    {
      info[j][symbol_index(line[j]) - 1]++;
      j++;
    }
    i++;
  }
}

/* create_cluster_files ***********************************************/
/*Donats dos indexs de clusters, crea un fitxer per cadascun, i escriu al
  fitxer les seqncies que el formen*/
void create_cluster_files(int i, int j)
{
  int k = 0;

  cluster_file_1 = fopen(cluster_filename_1, "w");
  cluster_file_2 = fopen(cluster_filename_2, "w");
  clusters_file = fopen(clusters_filename, "r");
  // Cluster i
  k = 0;
  while (k < clusters[i]->num_seqs)
  {
    pos_seq_clusters_file[seqs[i][k]] = add_sequence_cluster(cluster_file_1, seqs[i][k]);
    k++;
  }

  // Cluster j
  k = 0;
  while (k < clusters[j]->num_seqs)
  {
    pos_seq_clusters_file[seqs[j][k]] = add_sequence_cluster(cluster_file_2, seqs[j][k]);
    k++;
  }

  fclose(clusters_file);
  fclose(cluster_file_1);
  fclose(cluster_file_2);
}

/* find_nearest_clusters **************************************************/
/*Function that searches for the two closest sequences (or clusters). Return the
  result as a row (i) and a column (j). This function only loops
  half of the matrix, since it is symmetric.

  If it finds a zero it may mean:
    - is the similarity value
    - belong to the same cluster
*/

void find_nearest_clusters(int *i, int *j)
{
  int i_act = 0, j_act = 0, trobat = 0;

  // Let's look for the first one
  trobat = 0;
  for (i_act = 0; (i_act < num_seqs) && (trobat == 0); i_act++)
  {
    for (j_act = i_act + 1; (j_act < num_seqs) && (trobat == 0); j_act++)
    {
      if (matriu_puntuacions[i_act][j_act] != 0)
      {
        trobat = 1; // similarity value !=0
        *(i) = i_act;
        *(j) = j_act;
      }
      else
      {
        if (cluster_equivalent[j_act] != i_act)
        {
          trobat = 1; // do not belong to the same cluster
          *(i) = i_act;
          *(j) = j_act;
          matriu_puntuacions[i_act][j_act]; // I don't know what he's doing?
        }
      }
    }
  }

  i_act = *(i);
  j_act = *(j);
  while (i_act < num_seqs)
  {
    j_act = i_act + 1;
    while (j_act < num_seqs)
    {
      if (matriu_puntuacions[i_act][j_act] > matriu_puntuacions[*(i)][*(j)])
      {
        if (matriu_puntuacions[i_act][j_act] != 0)
        {
          *(i) = i_act;
          *(j) = j_act;
        }
        else
        {
          if (cluster_equivalent[j_act] != cluster_equivalent[i_act])
          {
            *(i) = i_act;
            *(j) = j_act;
          }
        }
      }
      j_act++;
    }
    i_act++;
  }
}

/* similarity_clusters *******************************************************/
// Calcula la similitud entre dos clusters. A la matriu
//"matriu" hi ha les similituds entre tot prefix dels dos clusters.
void similarity_clusters(float **matriu, long long_cluster1, long long_cluster2, float **info1, float **info2, long *puntuacio, int numseqs1, int numseqs2)
{
  int i = 0, j = 0;
  float p_diag = 0, p_esq = 0, p_amunt = 0, max = 0; // Puntuacions de les tres opcions possibles domplir una casella de la matriu
  char cami;

  // i->cluster1, j->cluster2

  matriu[0][0] = 0;
  matriu_cami[0][0] = 'd';
  for (i = 1; i <= long_cluster1; i++)
  {
    matriu[i][0] = matriu[i - 1][0] + calcul_simil_ini(info1[i - 1], numseqs2);
    matriu_cami[i][0] = 'a';
  }

  for (j = 1; j <= long_cluster2; j++)
  {
    matriu[0][j] = matriu[0][j - 1] + calcul_simil_ini(info2[j - 1], numseqs1);
    matriu_cami[i][0] = 'e';
  }

  for (i = 1; i <= long_cluster1; i++)
  {
    for (j = 1; j <= long_cluster2; j++)
    {
      p_diag = matriu[i - 1][j - 1] + calcul_simil(info1[i - 1], info2[j - 1]);
      p_esq = matriu[i][j - 1] + calcul_simil_ini(info2[j - 1], numseqs1);
      p_amunt = matriu[i - 1][j] + calcul_simil_ini(info1[i - 1], numseqs2);

      max = maxim_real(p_amunt, p_diag, p_esq, &cami);
      matriu[i][j] = max;
      matriu_cami[i][j] = cami;
    }
  }
  *puntuacio = max; // guarda la puntuacio del cluster
}

void mostrar_matriu_cars(char **m, int files, int cols)
{
  int i, j;

  for (i = 0; i < files; i++)
  {
    for (j = 0; j < cols; j++)
    {
      printf("%c ", m[i][j]);
    }
    printf("\n");
  }
}
/*******************************/

float calcul_simil_ini(float *freq, int n)
{
  // rep una frequencia de simbols i un enter simbolitzant gaps
  // es per inicialitzar la matriu del calcul de similituds entre clusters
  float sum;
  int i;

  sum = 0;
  // interclusters
  for (i = 0; i < num_symbols; i++)
    sum = sum + freq[i] * matpenal[i + 1][num_symbols] * n;
  // freq es des de 0 i marpenal des de 1
  sum = sum + intracluster(freq);
  return sum;
}
/*******************************/

float calcul_simil(float *freq1, float *freq2)
{
  // rep dues frequencies de simbols
  // es per inicialitzar la matriu del calcul de similituds entre clusters
  float sum;
  sum = 0;
  sum = sum + intracluster(freq1);
  sum = sum + intracluster(freq2);
  sum = sum + intercluster(freq1, freq2);
  return sum;
}
/****************************************/
float intercluster(float *freq1, float *freq2)
{
  float sum;
  int i, j;
  sum = 0;
  for (i = 0; i < num_symbols; i++)
    for (j = 0; j < num_symbols; j++)
      sum = sum + freq1[i] * freq2[j] * matpenal[i + 1][j + 1];
  return sum;
}

/****************************************/
float intracluster(float *freq)
{
  float sum;
  int i, j;
  sum = 0;
  // intracluster simbols diferents
  for (i = 0; i < num_symbols; i++)
    for (j = i + 1; j < num_symbols; j++)
      sum = sum + freq[i] * freq[j] * matpenal[i + 1][j + 1];
  // intracluster simbols iguals
  for (i = 0; i < num_symbols; i++)
    if (freq[i] > 1)
      sum = sum + freq[i] * (freq[i] - 1) * matpenal[i + 1][i + 1] / 2;
  return sum;
}

/*alignment_clusters ***********************************************/
/*Align clusters i, j*/
void alignment_clusters(float **matriu, int i, int j, float **info1, float **info2, char *line_orig, char *line_dest, char *resultat_cluster1, char *resultat_cluster2, long *puntuacio)
{
  int k = 0, mem_ok = 1;
  long res = 0;
  int len = 0, aux = 0;

  similarity_clusters(matriu, clusters[i]->long_seqs, clusters[j]->long_seqs, info1, info2,
                     puntuacio, clusters[i]->num_seqs, clusters[j]->num_seqs);
  // tenim la matriu de puntuacio i camins dels dos clusters ja calculada

  // printf("clusters: %d - %d\n", i, j);
  // mostrar_matriu(matriu, clusters[i]->long_seqs + 1, clusters[j]->long_seqs + 1);
  // mostrar_matriu_cars(matriu_cami, clusters[i]->long_seqs + 1, clusters[j]->long_seqs + 1);

  align_clusters(clusters[i]->long_seqs, clusters[j]->long_seqs, resultat_cluster1, resultat_cluster2, &len);
  rebuild_new_cluster(line_orig, line_dest, resultat_cluster1, resultat_cluster2, i, j, len);
  recalculate_matrix(i, j);
}

// rebuild_new_cluster ***********************************************/
// With a vector of ones and twos, and given the previous sequences, it gives the
// resulting cluster
void rebuild_new_cluster(char *line_orig, char *line_dest, char *resultat_cluster1, char *resultat_cluster2, int i, int j,
                             int len)
{
  // char *line_orig, *line_dest;
  int k = 0, l = 0, m = 0, n = 0;
  int seqs_new_cluster[300];

  clusters_file = fopen(clusters_filename, "a+");
  fseek(clusters_file, 0, SEEK_END); /*I DON'T KNOW IF IT IS NECESSARY ***************************************/
  /*cluster i*/
  k = 0;
  while (k < clusters[i]->num_seqs)
  {
    fseek(cluster_file_1, pos_seq_clusters_file[seqs[i][k]], SEEK_SET);
    fgets(line_orig, MAXLONGSEQ, cluster_file_1);
    l = 0; // pointer read array
    m = 0; // pointer writing array dest
    while (l < len)
    {
      if (resultat_cluster1[l] == '1') /*nucleotid*/
      {
        line_dest[l] = line_orig[m];
        line_dest[l + 1] = EOS;
        m++;
      }
      else
      { // gap
        line_dest[l] = '-';
        line_dest[l + 1] = EOS;
      }
      l++;
    }
    pos_seq[seqs[i][k]] = ftell(clusters_file);
    fprintf(clusters_file, "%s\n", line_dest);
    seqs_new_cluster[k] = seqs[i][k];
    k++;
  }
  n = clusters[i]->num_seqs;
  clusters[i]->num_seqs = clusters[i]->num_seqs + clusters[j]->num_seqs;
  clusters[i]->long_seqs = len;

  /*cluster j*/
  k = 0;
  while (k < clusters[j]->num_seqs)
  {
    fseek(cluster_file_2, pos_seq_clusters_file[seqs[j][k]], SEEK_SET);
    fgets(line_orig, MAXLONGSEQ, cluster_file_2);
    l = 0; // pointer read array
    m = 0; // pointer writing array dest
    while (l < len)
    {
      if (resultat_cluster2[l] == '1') // nucleotid
      {
        line_dest[l] = line_orig[m];
        line_dest[l + 1] = EOS;
        m++;
      }
      else
      { // gap
        line_dest[l] = '-';
        line_dest[l + 1] = EOS;
      }
      l++;
    }
    pos_seq[seqs[j][k]] = ftell(clusters_file);
    fprintf(clusters_file, "%s\n", line_dest);
    seq_cluster[seqs[j][k]] = i;
    seqs_new_cluster[n] = seqs[j][k];
    seqs_new_cluster[n + 1] = EOS;
    n++;
    k++;
  }
  fprintf(clusters_file, "%s\n", line_dest);

  seqs[i] = (int *)malloc(clusters[i]->num_seqs * sizeof(int) + 1);
  n = 0;
  while (n < clusters[i]->num_seqs)
  {
    seqs[i][n] = seqs_new_cluster[n];
    seqs[i][n + 1] = EOS;
    n++;
  }

  /* update the equivalent clusters of previous alignments */
  n = 0;
  while (n < num_seqs)
  {
    if ((cluster_equivalent[n] == j) || ((cluster_equivalent[n] == i) && (n != i))) // Just look at the "j". The "i" does not change
    {
      clusters[n] = clusters[i];
      seqs[n] = seqs[i];
      cluster_equivalent[n] = i;
    }
    n++;
  }
  fclose(clusters_file);
}
/**********************************************/
// from matriu_cami returns res1,res2 which are strips of {1,2}
//  where 1 indicates symbol and 2 gap

void align_clusters(int i, int j, char *res1, char *res2, int *len)
{
  if ((i == 0) && (j == 0))
  {
    *len = 0;
  }
  else
  {
    switch (matriu_cami[i][j])
    {
    case 'a':
      align_clusters(i - 1, j, res1, res2, len);
      (*len)++;
      res1[(*len) - 1] = '1'; // 1: there is a nucleotide in that position
      res2[(*len) - 1] = '2'; // 2: there is a gap in that position
      break;
    case 'd':
      align_clusters(i - 1, j - 1, res1, res2, len);
      (*len)++;
      res1[(*len) - 1] = '1';
      res2[(*len) - 1] = '1';
      break;
    case 'e':
      align_clusters(i, j - 1, res1, res2, len);
      (*len)++;
      res1[(*len) - 1] = '2';
      res2[(*len) - 1] = '1';
      break;
    default:
      break;
    }
  }
}

// add_sequence_cluster ******************************************/
long add_sequence_cluster(FILE *fitxer, int num_seq)
{
  long l = 0, pos_fitxer = 0;
  char seq[MAXLONGSEQ];

  if (pos_seq[num_seq] == -1)
  {
    // vol dir que la sequencia no s'ha escollit fins ara
    l = dona_longitud_seq(num_seq);
    clusters[num_seq]->long_seqs = l;

    carregar_sequencia(seq, num_seq);
    //    check_load_sequence(seq,clusters[num_seq]->long_seqs);
    pos_fitxer = ftell(fitxer);
    fprintf(fitxer, "%s", seq);
  }
  else
  {
    // it means it was included in another cluster
    fseek(clusters_file, pos_seq[num_seq], SEEK_SET);
    fgets(seq, MAXLONGSEQ, clusters_file);
    pos_fitxer = ftell(fitxer);
    fprintf(fitxer, "%s", seq);
  }
  return (pos_fitxer);
}

// recalculate_matrix ****************************************/
// Fix the matrix, knowing that we have aligned clusters "i" and "j"
void recalculate_matrix(int i, int j)
{
  int aux = 0, n = 0;

  aux = 0;
  while (aux < num_seqs)
  {
    if (i != aux)
    {
      if (cluster_equivalent[aux] == i)
      {
        matriu_puntuacions[i][aux] = 0;
        matriu_puntuacions[aux][i] = 0;
        matriu_puntuacions[j][aux] = 0;
        matriu_puntuacions[aux][j] = 0;
      }
      else
      {
        matriu_puntuacions[i][aux] = (matriu_puntuacions[i][aux] + matriu_puntuacions[j][aux]) / 2;
        matriu_puntuacions[aux][i] = matriu_puntuacions[i][aux];
        matriu_puntuacions[j][aux] = matriu_puntuacions[i][aux];
        matriu_puntuacions[aux][j] = matriu_puntuacions[i][aux];
      }
    }
    aux++;
  }

  aux = 0;
  while (aux < num_seqs)
  {
    if ((aux != i) && (aux != j))
    {
      if (cluster_equivalent[aux] == cluster_equivalent[i])
      {
        n = 0;
        while (n < num_seqs)
        {
          if ((aux != n) && (j != n) && (cluster_equivalent[aux] != cluster_equivalent[n]))
          {
            matriu_puntuacions[aux][n] = matriu_puntuacions[i][n];
            matriu_puntuacions[n][aux] = matriu_puntuacions[i][n];
          }
          else
          {
            if (cluster_equivalent[aux] == cluster_equivalent[i])
            {
              matriu_puntuacions[aux][n] = 0;
              matriu_puntuacions[n][aux] = 0;
            }
          }
          n++;
        }
      }
    }
    aux++;
  }
}
