OBJECTS = align.o ent_seq.o auxiliary.o pairs.o sort_seq.o multiple.o seqs_util.o checks.o

SOURCES = align.c ent_seq.c auxiliary.c pairs.c sort_seq.c multiple.c seqs_util.c checks.c

CFLAGS  = -c -ggdb

PROGNAME = alfm

alfm : $(OBJECTS)
	gcc -o $(PROGNAME) $(OBJECTS) -lm -lgd -lpng -ljpeg -lz -g -I/usr/include

align.o : align.c
	gcc $(CFLAGS) align.c

.c.o : $(SOURCES)
	gcc $(CFLAGS) $?

.PHONY: clean
clean :
	-rm $(PROGNAME) $(OBJECTS)