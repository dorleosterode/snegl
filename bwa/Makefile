BWA_DIR?=/usr/bin
CC=			gcc
CFLAGS=		-g -Wall -Wno-unused-function -O2 -I$(BWA_DIR)

smems: smems.o
		$(CC) $(CFLAGS) smems.o -o $@ -L$(BWA_DIR) -lbwa -lz

clean:
		rm -f gmon.out *.o a.out smems *~ *.a
