#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "bwt.h"
#include "bwa.h"
#include "kseq.h" // for the FASTA/Q parser
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

int main(int argc, char *argv[])
{
	bwaidx_t *idx;
	FILE *out;
	gzFile fp;
	kseq_t *ks;

	if (argc < 5) {
		fprintf(stderr, "Usage: smems <idx.base> <reads.fq|.fa> <output.txt> <min-length>\n");
		return 1;
	}

	char *filename = argv[2];
	char *outfile = argv[3];
	int min_len = atoi(argv[4]);
	if (min_len <= 0) {
	    fprintf(stderr, "Given min_length: %d. This is not a reasonable minimal length!\n", min_len);
	    fprintf(stderr, "Please enter a number between 1 and 2^32-1.\n");
	    exit(EXIT_FAILURE);
	}

	idx = bwa_idx_load(argv[1], BWA_IDX_ALL); // load the BWA index
	if (NULL == idx) {
	    fprintf(stderr, "Index load failed.\n");
	    exit(EXIT_FAILURE);
	}
	fp = gzopen(filename, "r");
	if (NULL == fp) {
	    fprintf(stderr, "Couldn't open %s\n", filename);
	    exit(EXIT_FAILURE);
	}
	out = fopen(outfile, "w");
	if (NULL == out) {
	    fprintf(stderr, "Couldn't open %s\n", outfile);
	    exit(EXIT_FAILURE);
	}

	ks = kseq_init(fp); // initialize the FASTA/Q parser

	while (kseq_read(ks) >= 0) { // read one sequence
	    int i, x = 0, k;
	    int start_width = 1;
	    int len = ks->seq.l;
	    uint8_t *seq = (uint8_t*) ks->seq.s;

	    for (i = 0; i < len; ++i) // convert to 2-bit encoding if we have not done so
		seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]];

	    bwtintv_v mem;
	    mem.m = 0;
	    mem.n = 0;
	    mem.a = 0;
	    //find all SMEMs
	    fprintf(out, ">%s %s\n", ks->name.s, ks->comment.s);
	    while (x < len) {
		if (seq[x] < 4) {
		    x = bwt_smem1(idx->bwt, len, seq, x, start_width, &mem, NULL);
		    for (i = 0; i < mem.n; ++i) {
			bwtintv_t p = mem.a[i];
			int slen = (uint32_t)p.info - (p.info>>32); // this is the length of the smem
			if (slen >= min_len) {
			    for (k = 0; k < p.x[2]; k++) {
				bwtint_t sfx_idx = p.x[0] + k;
				bwtint_t ref_pos = bwt_sa(idx->bwt, sfx_idx);
				fprintf(out, "%lu\t%u\t%u\n", ref_pos + 1, (uint32_t) (p.info >> 32) + 1, slen);
			    }
			}
		    }
		} else ++x;
	    }
	}

	kseq_destroy(ks);
	gzclose(fp);
	fclose(out);
	bwa_idx_destroy(idx);
	return 0;
}
