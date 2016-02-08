# snegl

SNEGL is a free software tool to compute maximal exact matches (MEMs)
between two dna or protein sequences. It can be used in a pipeline for
the computation of alignments based on a seed-and-extend
approach. SNEGL is based on an extension of the free software library
SDSL (see https://github.com/simongog/sdsl-lite/). The needed clone of
the SDSL can be found here
(https://github.com/dorleosterode/sdsl-lite/tree/dna_optimized_wt).

SNEGL can also compute maximal unique matches (MUMs) and super-maximal
exact matches (SMEMs), which can be used as seeds as well or as
representation of MEMs.

The directory bwa contains a test-programm for the computation of
SMEMs. It depends on the BWA library (see https://github.com/lh3/bwa).

Requirments and Dependencies
------------

snegl depends on a fork of the sdsl-lite library (see:
https://github.com/simongog/sdsl-lite, the fork in the right branch
can be found here:
https://github.com/dorleosterode/sdsl-lite/tree/dna_optimized_wt), so
you have to install it. You can follow the install instructions given
in the README. To use snegl you further have to fulfill the following
requirements:

* A modern, C++11 ready compiler such as `g++` version 4.7 or higher
  or `clang` version 3.2 or higher.

snegl is only tested on 64-bit Linux-systems, but you can try to use
it on different architectures.

The test-programm BWA-SMEM depends on the BWA library (see:
https://github.com/lh3/bwa). To use the test-programm you have to
install BWA using the install instructions in the README.

BWA-SMEM is also tested only on 64-bit Linux-systems.

Download and Compilation
------------

You can download snegl with the following command.

```sh
git clone https://github.com/dorleosterode/snegl.git
```

To compile snegl use:

```sh
make SDSL_DIR=/path/to/sdsl
```

The variable SDSL_DIR contains the absolute path to the directory,
where the sdsl-lite is installed. The default value for SDSL_DIR is
/usr/bin. If you have installed the sdsl-lite library there you can
just use:

```sh
make
```

If you need a static linked version of snegl use:

```sh
make static
```

You can remove the executable and all object-files with:

```sh
make clean
```

To compile BWA-SMEM use:
```sh
cd bwa
make BWA_DIR=/path/to/bwa
```

You can remove BWA-SMEM and all object-files with:
```sh
make clean
```

Getting started
----------

To compute the MEMs of minimal length 20 for a reference Sequenz
REF.fa and a query sequence QUERY.fa in the FASTA format use the
following command:

```sh
./snegl -l 20 REF.fa QUERY.fa
```

The results are written to a file called "output.sdsl.out". You can
change the filename with the parameter -o. So calling:

```sh
./snegl -l 20 -o new_output.txt REF.fa QUERY.fa
```

will write the results in the file "new_output.txt".

You can compute MUMs for the same setting as before with:

```sh
./snegl -l 20 --mum REF.fa QUERY.fa
```

and SMEMs with:

```sh
./snegl -l 20 --smem REF.fa QUERY.fa
```

respectivly. For other options see:

```sh
./snegl --help
```