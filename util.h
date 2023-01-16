/***************************************************************************
                                                            util.h  -  description
                                                               -------------------
                                      begin                : Tue Apr 24 2001
                                      copyright            : (C) 2001 by Roman Roset
                                      email                : roman.roset@menta.net

                                                           www.fib.upc.es
                                                       www.lsi.upc.es/~alggen
                                                         www.cepba.upc.es
                                   ***************************************************************************/

                                  /***************************************************************************
                                   *                                                                         *
                                   *  This program is property of U.P.C (Technical University of Catalonia)  *
                                   *  This program is made in CIRI (cepba-ibm ,European Center for           *
                                   *  Parallelism of Barcelona) like project for the Group ALGGEN            *
                                   *  (Algorithmics and Genetics Group).                                     *
                                   *                                                                         *
                                   **************************************************************************/

                                  #ifndef _DNASTICK_UTIL_H
                                  #define _DNASTICK_UTIL_H

                                  #include <time.h>
                                  #include <sys/time.h>
                                  #include <string.h>
                                  #include <stdlib.h>

                                  #define ASCII 0
                                  #define HTML 1

                                  #define RANDOM srand((unsigned int)time((time_t *)NULL));

                                  #define RANDOMINT(r) ((int)(rand() % (r)))

                                  #define min(a,b) (a<b?a:b)

                                  #define max(a,b) (a>b?a:b)

                                  #define GETTIME(time) \
                                          do { \
                                                  struct timeval _T_; \
                                                  gettimeofday(&_T_,NULL); \
                                                  time = _T_.tv_sec*1000 + _T_.tv_usec/1000; \
                                          } while (0)

                                  #define STRREV(seq) \
                                          do { \
                                                  int l = strlen(seq); \
                                                  char *_AUX_ = malloc(l + 1); \
                                                  int i = l; \
                                                  while (i) \
                                                              _AUX_[--i] = seq[l - i - 1]; \
                                                     _AUX_[l] = 0; \
                                                  free(seq); \
                                                  seq = _AUX_; \
                                          } while (0)

                                  #define CPCHARS(value, size) \
                                          (memcpy(malloc(size), value, size))

                                  #define CPVALUE(value, type) CPCHARS(value, sizeof (type))

                                  #define SUBSTR(res, str, from, to) \
                                          do { \
                                                  int  f = from ; \
                                                  int  t = to   ; \
                                                  res = calloc(t - f + 1, sizeof (char)); \
                                                  memcpy(res, str + f, t - f); \
                                          } while (0)

                                  #define FREE(elem) \
                                          if (elem) { \
                                                  free(elem); \
                                                  elem = NULL; \
                                          }
                                          
                                  #ifdef PARALLEL_OMP
                                          #define _INLINE_ _Inline
                                  #else
                                          #define _INLINE_
                                  #endif

                                  #define IO_NEW_LINE \
                                          (printf ("\n"))
                                          
                                  #define IO_PRINT_LINE(seq, nline) \
                                  do { \
                                          char *str = NULL; \
                                          SUBSTR(str , seq, nline * 80, (nline + 1) * 80);\
                                          printf("%s\n",str); \
                                  }while (0)        
                                          
                                  typedef enum {false=0, true} boolean;

                                  #endif
