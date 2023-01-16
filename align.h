/*************************************************************************/
                                  //Fitxer: align.h
                                  //Autor: Xavier Sol� Acha
                                  //Descripci�: Aquest fitxer contindr� totes les estructures de dades del 
                                  //programa, aix� com les diverses variables globals que necessitem.
                                  /*************************************************************************/
                                  #include <stdio.h>
                                  #include "const.h"
                                  FILE *fitxer_entrada; //Fitxer d'entrada de les seq��ncies
                                  FILE *fitxer_temporal; //Fitxer temporal que emmagatzema les seq��ncies
                                  FILE *fitxer_alfabet;//Fitxer que emmagatzema l'alfabet en el formata adequat

                                  char nom_fitxer_temporal[13]; //Nom del fitxer temporal

                                  int numsimb;   //nombre de simbols de l'alfabet
                                  char alfabet[MAXLONGALF];//alfabet en posicions 1,2,..,numsimb
                                  int num_seqs;   //Nombre de seq��ncies que hi ha al fitxer

                                  float matpenal[MAXLONGALF][MAXLONGALF];  //PUNTUACIO match,mismatch i gap...
                                  float **matriu_puntuacions;  //Matriu que contindr� la puntuaci� de l'alineament �ptim entre totes les seq��ncies
                                  char **matriu_cami; //Cont� per on s'ha omplert cada posicio de la matriu ('e': esquerra, 'a':amunt, 'd': diagonal)

                                  int args(int argc, char** argv);
                                  void leer_alfabeto(FILE *fd);                     
