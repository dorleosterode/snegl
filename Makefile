CC= g++
CFLAGS=	-ggdb -std=c++11 -Wall -O3 -DNDEBUG
LDFLAGS= -static
PROG= snegl
SDSL_DIR?= /usr/bin
INCLUDE= $(SDSL_DIR)/include
LIBSINCLUDE= $(SDSL_DIR)/lib
LIBS= -lsdsl -ldivsufsort -ldivsufsort64 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
SUBDIRS=	.

.SUFFIXES:.c .o .cpp

.c.o:
		$(CC) -c $(CFLAGS) -I$(INCLUDE) $< -o $@
.cpp.o:
		$(CC) -c $(CFLAGS) -I$(INCLUDE) $< -o $@

all:$(PROG)

snegl:process_text.o utility.o dna_algorithms.o process_dna.o mem.o
		$(CC) $(CFLAGS) process_text.o utility.o dna_algorithms.o process_dna.o mem.o -o $@ -L$(LIBSINCLUDE) $(LIBS)

static:process_text.o utility.o dna_algorithms.o process_dna.o mem.o
		$(CC) $(CFLAGS) $(LDFLAGS) process_text.o utility.o dna_algorithms.o process_dna.o mem.o -o $@ -L$(LIBSINCLUDE) $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.cpp )

# DO NOT DELETE THIS LINE -- make depend depends on it.

process_text.o: utility.hpp process_text.hpp sdsl_types.hpp
dna_algorithms.o: sdsl_types.hpp dna_algorithms.hpp
process_dna.o: process_dna.hpp sdsl_types.hpp dna_algorithms.hpp utility.hpp
mem.o: process_text.hpp process_dna.hpp sdsl_types.hpp
utility.o: utility.hpp sdsl_types.hpp
