/*************************************************************************/
// Fitxer: align.c
// Autor: Xavier Sol Acha
// Descripci: Aquest fitxer contindr el programa principal de l'aplicaci
/*************************************************************************/

#include "align.h"
#include "auxiliar.h"
#include "ent_seq.h"
#include "parells.h"
#include "const.h"
#include "multiple.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <ctype.h>
#include "comprobacions.h"

FILE *fitxer_entrada;  // Fitxer d'entrada de les seqncies
FILE *fitxer_temporal; // Fitxer temporal que emmagatzema les seqncies
FILE *fitxer_alfabet;  // Fitxer que emmagatzema l'alfabet en el formata adequat

char nom_fitxer_temporal[13]; // Nom del fitxer temporal

int numsimb;              // nombre de simbols de l'alfabet
char alfabet[MAXLONGALF]; // alfabet en posicions 1,2,..,numsimb
int num_seqs;             // Nombre de seqncies que hi ha al fitxer

float matpenal[MAXLONGALF][MAXLONGALF]; // PUNTUACIO match,mismatch i gap...
float **matriu_puntuacions;             // Matriu que contindr la puntuaci de l'alineament ptim entre totes les seqncies
char **matriu_cami;                     // Cont per on s'ha omplert cada posicio de la matriu ('e': esquerra, 'a':amunt, 'd': diagonal)

int args(int argc, char **argv)
{
  // comprueba los argumentos
  int i;
  if (argc < 3)
  {
    fprintf(stderr, "Sintaxis: %s fitxeralfabet fitxersequenciesformatfasta \n", argv[0]);
    printf("El format del fitxer alfabet ha de ser:\n");
    printf("    Numero de simbols\n");
    printf("    Simbols de l'alfabet separats per un caracter blanc\n");
    printf("    Penalitzacio del gap(nomes n'hi ha d'un tipus)\n");
    printf("i les seguents linies 1,2,...,i, contenen 1 numero, 2 numeros,...\n");
    printf("que son les puntuacions del simbol i el gap amb els simbols anteriors");
    printf("EXEMPLE de 6 simbols (inclos el gap)\n");
    printf("6\n");
    printf("o p s n c -\n");
    printf(" 2     Puntuacio primer simbol amb ell matex\n");
    printf("-1  2    Puntuacio segon simbol amb primer i ell mateix\n");
    printf("-2 -2  1   ...\n");
    printf("-2 -2  0  1\n");
    printf("-2 -2 -1 -1 1\n");
    printf("-2 -2  0  0 0 0\n");
    return (1);
  }
  else
    return (0);
}

// funcio que llegeix l'alfabet,la matriu de puntuacio i el gap
void leer_alfabeto(FILE *fd)
{
  int i, j;
  fscanf(fd, "%d\n", &numsimb);

  for (i = 1; i <= numsimb; i++)
  {
    fscanf(fd, "%c ", &alfabet[i]);
  }
  for (i = 1; i < numsimb; i++)
  {
    alfabet[i] = toupper(alfabet[i]);
  }
  for (i = 1; i <= numsimb; i++)
  {
    for (j = 1; j <= i; j++)
    {
      fscanf(fd, "%f", &matpenal[i][j]);
      matpenal[j][i] = matpenal[i][j];
    }
  }
}

int main(int argc, char *argv[])
{
  int i = 0, j, k, l, estat;
  long long_seq, long_seq2;
  char *seq, *seq2;
  float **matriu, res;
  char pref[6];
  int num_pref;

  if (args(argc, argv) == 1)
    exit(-1);
  if ((fitxer_alfabet = fopen(argv[1], "r")) == NULL)
  {
    printf("Error en lectura del fichero alfabeto\n");
    exit(1);
  }
  leer_alfabeto(fitxer_alfabet);
  comprovarlecturaalfabet();
  // inicialitza les variables globals alfabet,matpenal,numsimb
  if ((fitxer_entrada = fopen(argv[2], "r")) == NULL)
  {
    printf("Error en lectura del fichero de secuencias\n");
    exit(1);
  }

  // Esborrem els fitxers temporals. Tanquem els que no necessitem en aquest moment
  num_pref = generar_prefix_fitxers(argv[1]);
  omplir_string_prefix(pref, num_pref);
  sprintf(nom_fitxer_temporal, "%s1.tmp", pref);
  sprintf(nom_fitxer_clusters, "%s2.tmp", pref);
  sprintf(nom_fitxer_cluster_1, "%s3.tmp", pref);
  sprintf(nom_fitxer_cluster_2, "%s4.tmp", pref);
  fitxer_temporal = fopen(nom_fitxer_cluster_1, "w");
  fclose(fitxer_temporal);
  fitxer_temporal = fopen(nom_fitxer_cluster_2, "w");
  fclose(fitxer_temporal);
  fitxer_temporal = fopen(nom_fitxer_clusters, "w");
  fclose(fitxer_temporal);
  fitxer_temporal = fopen(nom_fitxer_temporal, "w");
  if (llegir_sequencies_fitxer() == -1)
  // te les sequencies en el fitxer temporal en un format especial
  // numeroseq,nom,longitud,sequencia
  {
    printf("Too many (>%d) or too few (<2) input sequences\n", MAXSEQ);
    remove(nom_fitxer_clusters);
    remove(nom_fitxer_cluster_1);
    remove(nom_fitxer_cluster_2);
    remove(nom_fitxer_temporal);
    exit(-1);
  }
  else // nombre de sequencies correcte
  {

    fclose(fitxer_temporal);
    //       comprovar_fitxer_temporal();
    fitxer_temporal = fopen(nom_fitxer_temporal, "r");

    estat = 0;
    matriu = (float **)malloc(MAXLONGALIN * sizeof(float *));
    if (matriu == NULL)
    {
      estat = -1;
      printf("Out of memory\n");
      remove(nom_fitxer_clusters);
      remove(nom_fitxer_cluster_1);
      remove(nom_fitxer_cluster_2);
      remove(nom_fitxer_temporal);
      exit(-1);
    }
    // hi ha prou memoria per a 2000 elements
    i = 0;
    estat = 0;
    while ((i < MAXLONGALIN) && (estat == 0))
    {
      matriu[i] = (float *)malloc(MAXLONGALIN * sizeof(float));
      if (matriu[i] == NULL)
      {
        estat = -1;
        printf("Out of memory\n");
        remove(nom_fitxer_clusters);
        remove(nom_fitxer_cluster_1);
        remove(nom_fitxer_cluster_2);
        remove(nom_fitxer_temporal);
        exit(-1);
      }
      i++;
    }
    // hi ha prou memoria per a una matriu de 2000x2000

    estat = 0;
    matriu_cami = (char **)malloc(MAXLONGALIN * sizeof(char *));
    if (matriu_cami == NULL)
    {
      estat = -1;
      printf("Out of memory\n");
      remove(nom_fitxer_clusters);
      remove(nom_fitxer_cluster_1);
      remove(nom_fitxer_cluster_2);
      remove(nom_fitxer_temporal);
      exit(-1);
    }
    i = 0;
    estat = 0;
    while ((i < MAXLONGALIN) && (estat == 0))
    {
      matriu_cami[i] = (char *)malloc(MAXLONGALIN * sizeof(char));
      if (matriu_cami[i] == NULL)
      {
        estat = -1;
        printf("Out of memory\n");
        remove(nom_fitxer_clusters);
        remove(nom_fitxer_cluster_1);
        remove(nom_fitxer_cluster_2);
        remove(nom_fitxer_temporal);
        exit(-1);
      }
      i++;
    }

    matriu_puntuacions = (float **)malloc(num_seqs * sizeof(float *));
    if (matriu_puntuacions == NULL)
    {
      estat = -1;
    }
    i = 0;
    estat = 0;
    while ((i < num_seqs) && (estat == 0))
    {
      matriu_puntuacions[i] = (float *)malloc(num_seqs * sizeof(float));
      if (matriu_puntuacions[i] == NULL)
      {
        estat = -1;
        printf("Out of memory\n");
        remove(nom_fitxer_clusters);
        remove(nom_fitxer_cluster_1);
        remove(nom_fitxer_cluster_2);
        remove(nom_fitxer_temporal);
        exit(-1);
      }
      i++;
    }

    j = 0;
    seq = (char *)malloc(MAXLONGSEQ * sizeof(char));
    seq2 = (char *)malloc(MAXLONGSEQ * sizeof(char));

    // hi ha prou memoria per a totes les variables
    // printf("\n 1");
    while (j < (num_seqs - 1))
    {
      long_seq = dona_longitud_seq(j);
      carregar_sequencia(seq, j);
      // comprovar_carregar_sequencia(seq,long_seq);
      k = j + 1;
      printf("\n");
      while (k < num_seqs)
      {
        printf(".");
        long_seq2 = dona_longitud_seq(k);
        carregar_sequencia(seq2, k);
        // comprovar_carregar_sequencia(seq2,long_seq2);
        res = similitud(matriu, seq, seq2, long_seq, long_seq2);
        //      comprova_matriu_cami(matriu,long_seq,long_seq2);
        //      comprovar_similitud(seq, seq2, long_seq, long_seq2);
        matriu_puntuacions[j][k] = res;
        matriu_puntuacions[k][j] = res;
        k++;
      }
      j++;
    }
    comprovar_matriu_similaritats(matriu_puntuacions);
    free(seq);
    free(seq2);
    // printf("\n 2");
    alineament_multiple(matriu, atoi(argv[2]));
    // segon argument es el format de sortida
    // printf("\n 3");
    i = 0;
    while (i < num_seqs)
    {
      free(matriu_puntuacions[i]);
      i++;
    }
    free(matriu_puntuacions);

    i = 0;
    while (i < MAXLONGALIN)
    {
      free(matriu[i]);
      i++;
    }
    free(matriu);

    i = 0;
    while (i < MAXLONGALIN)
    {
      free(matriu_cami[i]);
      i++;
    }
    free(matriu_cami);
    fclose(fitxer_temporal);
    remove(nom_fitxer_clusters);
    remove(nom_fitxer_temporal);
    remove(nom_fitxer_cluster_1);
    remove(nom_fitxer_cluster_2);
  }
}
