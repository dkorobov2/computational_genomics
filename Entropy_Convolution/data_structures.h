#ifndef _DATA_STRUCTURES_H
#define _DATA_STRUCTURES_H

//Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <bitset>
#include <random>
#include <bitset>
#include "RTree.h"
#include "params.h"

using namespace std;

typedef struct ncount_t {
	uint32_t num_a;
	uint32_t num_c;
	uint32_t num_g;
	uint32_t num_t;
} ncount_t;

typedef struct rtree_node {
	int start;
	int finish;
} rtree_node;

typedef RTree<int, int, 1, float> SRTree;

ncount_t * createSampleCounts(int length);

uint32_t * createSampleReference(int length);

/* Creates a mock array with starting positions. First element is the total number of short reads */
uint32_t * createStarting(int short_reads);

/* Creates a mock array with short reads. First element is the total number of short reads */
uint32_t * createNucs(int short_reads);

/* Creates a mock array of cigar strings. First element is the total numer of short reads */
uint32_t * createCigars(int len);

/* CPU version of entropy calculation function */
void calculateEntropy(float * output, ncount_t * counts, uint32_t * reference);

/* Returns the two bits corresponding to the nucleotide in the reference at that position */
uint32_t nucleotide(uint32_t * reference, int position);

/* Printing functions */
void printEntropy(float * entropy, int num_elem);

void printConvolution(float * convolution, int num_elem);

void printThreshold(unsigned char * threshold, int num_elem);

void printCounts(ncount_t * counts, int num_elem);

void print_elements(vector<uint32_t>, int n);

/* Convolution function on the CPU for testing */
void convolve(float Signal[/* SignalLen */], size_t SignalLen,
	float Kernel[/* KernelLen */], size_t KernelLen,
	float Result[/* SignalLen + KernelLen - 1 */]);

// Call back function for RTree search
bool MySearchCallback(int id, void* arg);

/* Used after the threshold kernel to create input used in the de bruijn assembly */
std::vector<std::uint32_t> createDBInput(uint32_t * const nuc, vector<unsigned char> threshold, vector<uint32_t> const &sr_indexes, SRTree &sr_tree);

/* Create R-tree with start and end positions of each short read */
void createRTreeAndFrequencies(uint32_t * const nuc, uint32_t * const starting, uint32_t * const cigar, SRTree &sr_tree, vector<ncount_t> &freq, vector<uint32_t> &sr_indexes);

/* Creates the array of nucleotide frequencies at each location */
std::vector<ncount_t> createNucFrequencies(uint32_t * const nuc, uint32_t * const starting, uint32_t * const cigar);

/* Computes the frequencies for an individual short read */
void computeFrequencies(uint32_t * const nuc, uint32_t * const cigar, vector<ncount_t> &freq, int start_index, int read_index);

/* Helper function used for summing the insertion and deletions to determine the correct end length */
int computeEndIndex(uint32_t * const nuc, uint32_t * const cigar, int start_index, int read_index);

/* Combines computeEndIndex and computeFrequencies into one function to reduce runtime */
int computeEndIndexAndFrequencies(uint32_t * const nuc, uint32_t * const cigar, vector<ncount_t> &freq, int start_index, int read_index);

/* Traverse the threshold array and determine if the short reads fall in a critical section */
int checkRead(uint8_t * const threshold, int start_index, int end_index);

/* Array used for calculating the extra elements needed to pass into the convolution kernel */
int * createExtraArray(int k, int grain_size);
#endif