/*
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cinttypes>
#include "mapping.h"

void addhash(uint32_t kmer, kmer_t **table, uint32_t read, uint32_t k_num) {
	kmer_t *temp = new kmer_t{kmer, read, k_num, NULL};
	uint8_t hash = kmer & 0b111111;
	if (table[hash] == NULL) {
		table[hash] = temp;
	}
	else {
		kmer_t * col = table[hash];
		temp->next = col;
		table[hash] = temp;
	}
}
*/