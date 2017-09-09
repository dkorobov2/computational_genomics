#ifndef _READS_H
#define _READS_H

#define READ_SIZE 36 //number of bytes in a single read
#define LENGTH 294912 //There are 32 * 1024 reads in this sequence (1152 * 1024 bytes).
#define NUM_NUCL 144 // # of "nucleotides" (including header) in a read

typedef struct kmer_t {
	uint32_t kmer;
	uint32_t k_num;
} kmer_t;

uint32_t * createSampleReads(int len);

uint32_t * createSmallSample();

void print_two_reads(uint32_t * sample);

void print_elements(uint32_t * sample, int n);

void print_freq_table(kmer_t * frequency_table, int n);

#endif