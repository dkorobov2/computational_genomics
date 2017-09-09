#ifndef _PARAMS_H
#define _PARAMS_H

#define ENTROPY "entropy.txt"
#define CONVOLUTION "convolution.txt"
#define THRESHOLD "threshold.txt"

#define NUM_WORKGROUPS 1024 // the number of workgroups to be processed at a time

#define LOCAL 512 //local size (max on my computer)
#define GRAIN_SIZE 1 //number of elements per thread (2), found this to be most efficient through testing
#define NUM_ITERATIONS 2 //arbitrary for now. number of times we iterate through the kernels
#define NUM_NUCL LOCAL * GRAIN_SIZE * NUM_WORKGROUPS // number of elements per iteration
#define NUM_NUCL_TOTAL LOCAL * GRAIN_SIZE * NUM_WORKGROUPS * NUM_ITERATIONS //total number of elements accounting for iterations

#define FILTER_LENGTH 5
#define THRESHOLD_VAL 6 //arbitrary for now

#define NUM_SR_ELEMS 9 // number of 4-byte elements in a short read for my samples
#define NUM_READS 10 // something for testing
#define NUC_PER_ELEM 16 // number of nucleotides per 4 byte element

#endif