OBJECTS = align.o ent_seq.o auxiliar.o parells.o sort_seq.o multiple.o seqs_util.o comprobacions.o

SOURCES = align.c ent_seq.c auxiliar.c parells.c sort_seq.c multiple.c seqs_util.c comprobacions.c

CFLAGS  = -c -ggdb

PROGNAME = alfm

alfm : $(OBJECTS)
	gcc -o $(PROGNAME) $(OBJECTS) -lm -lgd -lpng -ljpeg -lz -g -I/usr/include

align.o : align.c
	gcc $(CFLAGS) align.c

.c.o : $(SOURCES)
	gcc $(CFLAGS) $?
