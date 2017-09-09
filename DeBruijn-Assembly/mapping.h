#ifndef _MAPPING_H
#define _MAPPING_H

//Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <math.h>
#include "reads.h"
#include <chrono>
using namespace std;

#define TOTAL_BITS 288 //number of bits
#define READ_BITS 256 //number of bits
#define LOCAL_BYTES LOCAL_SIZE*4 //number of bytes in a work group (1440)

#define NUM_NUCL 144 // # of "nucleotides" (including header) in a read

typedef pair<uint32_t, uint32_t> h_pair;

void createStdMap(uint32_t *read_seq, int k, multimap<uint64_t, h_pair> * hashmap, int seq_len);

void createMap(uint32_t *read_seq, int k, uint32_t *table, int seq_len);

unsigned int nextPowerOf2(unsigned int n);

void countFrequency(uint32_t * table_in, kmer_t * table_out, int size);

#endif