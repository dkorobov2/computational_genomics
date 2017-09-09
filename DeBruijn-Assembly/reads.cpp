//Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "reads.h"

using namespace std;

/*  Creates a sample sequence of reads of size READ_SIZE * NUM_READS */
uint32_t * createSampleReads(int len) {
	if (len % 9 != 0) {
		return NULL;
	}
	uint32_t test1 = 0xFFFFFFFF;

	uint32_t test2 = 0xFFFFFFFF;

	uint32_t * read_seq = new uint32_t[len];
	for (int i = 0; i < len; i++) {
		if (i % 9 == 0) {
			read_seq[i] = 128;
		}
		else {
			if (i % 2)
				read_seq[i] = test1;
			else
				read_seq[i] = test2;
		}
	}
	return read_seq;
}

/*  Creates small sample sequence of reads (2 reads) */
uint32_t * createSmallSample() {
	uint32_t test1 = 0xAAAAAAAA;

	uint32_t test2 = 0x55555555;

	uint32_t * read_seq = new uint32_t[READ_SIZE * 2];
	for (int i = 0; i < READ_SIZE * 2; i++) {
		if (i % 9 == 0) {
			read_seq[i] = 128;
		}
		else {
			if (i % 2)
				read_seq[i] = test1;
			else
				read_seq[i] = test2;
		}
	}
	return read_seq;
}

/* prints first 2 reads in a sequence */
void print_two_reads(uint32_t * sample) {
	for (int i = 0; i < 18; i++) {
		cout << hex << sample[i] << endl;
	}
}

/* prints first n integer elements in hex */
void print_elements(uint32_t * sample, int n) {
	for (int i = 0; i < n; i++) {
		cout << hex << sample[i] << endl;
	}
}

void print_freq_table(kmer_t * frequency_table, int n) {
	for (int i = 0; i < 100; i++) {
		printf("Kmer: %x           Frequency: %d\n", frequency_table[i].kmer, frequency_table[i].k_num);
	}
}