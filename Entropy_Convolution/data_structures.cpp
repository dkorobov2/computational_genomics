#include "data_structures.h"

template< size_t size>
typename std::bitset<size> random_bitset(double p = 0.5) {

	typename std::bitset<size> bits;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::bernoulli_distribution d(p);

	for (int n = 0; n < size; ++n) {
		bits[n] = d(gen);
	}

	return bits;
}

ncount_t * createSampleCounts(int length)
{
	ncount_t *counts = new ncount_t[length];
	for (int i = 0; i < length; i++) {
		counts[i].num_a = rand() % 25 + 1;
		counts[i].num_c = rand() % 25 + 1;
		counts[i].num_g = rand() % 25 + 1;
		counts[i].num_t = rand() % 25 + 1;
	}
	return counts;
}

// not being used
uint32_t * createSampleReference(int length)
{
	uint32_t *reference = new uint32_t[length];
	for (int i = 0; i < length; i++) {
		reference[i] = random_bitset<32>().to_ulong();
		//bitset<32> x(reference[i]);
		//cout << x << endl;
	}
	return reference;
}

// This is the starting array, length represents the number of short reads
uint32_t * createStarting(int short_reads)
{
	uint32_t *starting = new uint32_t[short_reads + 1];
	starting[0] = short_reads;
	for (int i = 1; i < short_reads + 1; i++) {
		// starting indexes are from 0 to 5000
		starting[i] = i - 1;
		//bitset<32> x(reference[i]);
		//cout << x << endl;
	}
	return starting;
}

// Creates array where first element is the number of reads, followed by regular short reads
// len is the number of short reads, array size will be 1 + 9 * len
uint32_t * createNucs(int short_reads) {
	uint32_t * nuc_seq = new uint32_t[1 + NUM_SR_ELEMS * short_reads];
	uint32_t * starting_index = &nuc_seq[1];
	nuc_seq[0] = short_reads;
	for (int i = 0; i < NUM_SR_ELEMS * short_reads; i++) {
		if (i % 9 == 0) {
			//read_seq[i] = rand() % 128;
			starting_index[i] = 128;
		}
		else {
			//read_seq[i] = rand() % UINT_MAX;
			starting_index[i] = 0xbbbbbbbb;
		}
	}
	return nuc_seq;
}

// Creates array where first element is the number of reads, followed by cigar strings
// length is the total number of uint32_t elements
// MATCH : 00
// DELETION : 01
// INSERTION : 10
// MISMATCH : 11
// len is the number of short reads, so the array size will be 1 + 9 * len
uint32_t * createCigars(int len) {
	uint32_t * cigar_seq = new uint32_t[1 + NUM_SR_ELEMS * len];
	uint32_t * starting_index = &cigar_seq[1];
	cigar_seq[0] = len;
	for (int i = 0; i < NUM_SR_ELEMS * len; i++) {
		if (i % 9 == 0) {
			//read_seq[i] = rand() % 128;
			starting_index[i] = 128;
		}
		else {
			starting_index[i] = 0;
			//starting_index[i] = 0x55555555;
			//starting_index[i] = random_bitset<32>().to_ulong();
		}
	}
	return cigar_seq;
}

void calculateEntropy(float * output, ncount_t * counts, uint32_t * reference)
{
	int length = sizeof(counts) / sizeof(counts[0]);
	int sum;
	float pa, pc, pt, pg;
	for (int i = 0; i < length; i++) {
		sum = counts[i].num_a + counts[i].num_c + counts[i].num_g + counts[i].num_t;
		switch (nucleotide(reference, i)) {
		case 0:
			pa = (float)counts[i].num_a / sum;
			output[i] = -sum * pa * log2(pa);
			break;
		case 1:
			pc = (float)counts[i].num_c / sum;
			output[i] = -sum * pc * log2(pc);
			break;
		case 2:
			pg = (float)counts[i].num_g / sum;
			output[i] = -sum * pg * log2(pg);
			break;
		case 3:
			pt = (float)counts[i].num_t / sum;
			output[i] = -sum * pt * log2(pt);
			break;
		}
	}
}

uint32_t nucleotide(uint32_t * reference, int position)
{
	int index = position / 16;
	int shift = (15 - (position % 16))*2;
	int nuc = (reference[index] >> shift) & 0x3;
	//printf("reference: %x shift: %d nuc: %d\n", reference[index], shift, nuc);
	return nuc;
}

void printCounts(ncount_t * counts, int num_elem)
{
	cout << "**********COUNTS**********" << endl;
	for (int i = 0; i < num_elem; i++) {
		cout << "A: " << counts[i].num_a << " C: " << counts[i].num_c << " G: " << counts[i].num_g << " T: " << counts[i].num_t << endl;
	}
}

void printEntropy(float * entropy, int num_elem)
{
	int length = sizeof(entropy) / sizeof(entropy[0]);
	cout << "**********ENTROPY**********" << endl;
	for (int i = 0; i < num_elem; i++) {
		cout << "entropy: " << entropy[i] << endl;
	}
}

void printConvolution(float * convolution, int num_elem)
{
	int length = sizeof(convolution) / sizeof(convolution[0]);
	cout << "**********CONVOLUTION**********" << endl;
	for (int i = 0; i < num_elem; i++) {
		cout << "convolution: " << convolution[i] << endl;
	}
}

void printThreshold(unsigned char * threshold, int num_elem)
{
	int length = sizeof(threshold) / sizeof(threshold[0]);
	cout << "**********THRESHOLD**********" << endl;
	for (int i = 0; i < num_elem; i++) {
		cout << "threshold: " << (int)threshold[i] << endl;
	}
}

/* prints first n integer elements in hex */
void print_elements(vector<uint32_t> sample, int n) {
	for (int i = 0; i < n; i++) {
		cout << hex << sample[i] << endl;
	}
}

void convolve(float Signal[/* SignalLen */], size_t SignalLen,
	float Kernel[/* KernelLen */], size_t KernelLen,
	float Result[/* SignalLen + KernelLen - 1 */])
{
	size_t n;

	for (n = 0; n < SignalLen + KernelLen - 1; n++)
	{
		size_t kmin, kmax, k;

		Result[n] = 0;

		kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
		kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

		for (k = kmin; k <= kmax; k++)
		{
			Result[n] += Signal[k] * Kernel[n - k];
		}
	}
}

// Call back function for RTree search
bool MySearchCallback(int id, void* arg)
{
	cout << "Hit short read number " << id << "\n";
	int * temp = (int *)arg;
	// initialize the min and max
	if (temp[0] == 0 && temp[1] == 0) {
		temp[0] = id;
		temp[1] = id;
	}
	//update the min
	else if (id < temp[0]) {
		temp[0] = id;
	}
	// update the max
	else if (id > temp[1]) {
		temp[1] = id;
	}
	return true; // keep going
}

std::vector<std::uint32_t> createDBInput(uint32_t * const nuc, vector<unsigned char> threshold, vector<uint32_t> const &sr_indexes, SRTree &sr_tree)
{
	vector<uint32_t> result;
	// the start and end indexes in the threshold array
	int start_index, end_index;
	uint32_t minmax[2];
	int num_reads = nuc[0];
	int min_read, max_read, num_elements, db_index = 0;
	uint8_t bool_start = 0;
	cout << "THRESHOLD SIZE: " << threshold.size() << endl;
	for (int i = 0; i < (int)threshold.size(); i++) {
		// if a string of ones is encountered
		if (threshold[i] == 1 && bool_start == 0) {
			start_index = i;
			bool_start = 1;
		}
		// middle element in string of 1's
		else if (threshold[i] == 1 && bool_start == 1) {
			// keep iterating
		}
		// if a 0 is found, end index is i-1
		else if (threshold[i] == 0 && bool_start == 1) {
			//reset the values of minmax
			minmax[0] = 0;
			minmax[1] = 0;
			// set the value of end_index
			end_index = i - 1;
			// reset bool_start 
			bool_start = 0;
			// find min_read and max_read using RTree search
			sr_tree.Search(&start_index, &end_index, MySearchCallback, minmax);
			min_read = minmax[0];
			max_read = minmax[1];
			if (max_read + 1 < num_reads) {
				// calculate the number of bytes to copy for db assembly input
				num_elements = sr_indexes[max_read + 1] - sr_indexes[min_read];
			}
			else if (max_read + 1 == num_reads) {
				num_elements = sizeof(&nuc[sr_indexes[min_read]]);
			}
			result.resize(db_index + num_elements + 1);
			memcpy(&result[db_index], &nuc[sr_indexes[min_read]], num_elements * sizeof(uint32_t));
			db_index += num_elements;
		}
	}
	return result;
}

void createRTreeAndFrequencies(uint32_t * const nuc, uint32_t * const starting, uint32_t * const cigar, SRTree &sr_tree, vector<ncount_t> &freq, vector<uint32_t> &sr_indexes)
{
	int num_reads = nuc[0];
	sr_indexes[0] = num_reads;
	cout << "num_reads: " << num_reads << endl;
	int read_index = 1;
	int start_index, end_index, sr_length;
	for (int i = 0; i < num_reads; i++) {
		start_index = starting[i + 1];
		sr_length = nuc[read_index];
		end_index = computeEndIndexAndFrequencies(nuc, cigar, freq, start_index, read_index);
		//cout << "start index: " << start_index << " end_index: " << end_index << "i: " << i << " read_index: " << read_index << endl;
		sr_tree.Insert(&start_index, &end_index, i);
		sr_indexes[i + 1] = read_index;
		// this is the length of the short read (assuming 32 bit elements)
		read_index += 1 + (int)ceil((end_index - start_index) / 16);
	}
}

/* This is no longer being used so we don't have to traverse the reads twice, once here and once for the RTree */
std::vector<ncount_t> createNucFrequencies(uint32_t * const nuc, uint32_t * const starting, uint32_t * const cigar)
{
	vector<ncount_t> result;
	int num_reads = nuc[0];
	int start_index;
	for (int i = 0; i < num_reads; i++) {
		start_index = starting[i + 1];
		computeFrequencies(nuc, cigar, result, start_index, i);
	}
	return result;
}

// NO LONGER BEING USED
void computeFrequencies(uint32_t * const nuc, uint32_t * const cigar, vector<ncount_t> &freq, int start_index, int read_index)
{
	uint32_t n, c, position;
	int sr_length = nuc[read_index];
	uint32_t * sr_start = &nuc[1 + read_index];
	uint32_t * cigar_start = &nuc[1 + read_index];
	position = 128 - sr_length;
	//cout << "vector size: " << freq.size() << endl;
	for (int i = 128 - sr_length; i < sr_length; i++) {
		n = nucleotide(sr_start, i);
		c = nucleotide(cigar_start, i);
		cout << "start index: " << start_index << " position " << position << endl;
		// if match or mismatch, increment the corresponding nucleotide count
		if (c == 0 || c == 3) {
			switch (n) {
			case 0:
				freq[start_index + position].num_a++;
				break;
			case 1:
				freq[start_index + position].num_c++;
				break;
			case 2:
				freq[start_index + position].num_g++;
				break;
			case 3:
				freq[start_index + position].num_t++;
				break;
			}
			position++;
		}
		// if insertion, don't increment the position and do nothing
		else if (c == 2) {
			/* do nothing */
		}
		// if deletion, just increment the position
		else if (c == 1) {
			position++;
		}
	}
}

// NO LONGER BEING USED
int computeEndIndex(uint32_t * const nuc, uint32_t * const cigar, int start_index, int read_index)
{
	uint32_t n;
	int num_deletions = 0;
	int sr_length = cigar[read_index];
	uint32_t * sr_start = &cigar[1 + read_index];
	for (int i = 128 - sr_length; i < sr_length; i++) {
		n = nucleotide(sr_start, i);
		if (n == 1)
			num_deletions++;
	}
	if (num_deletions > 0) {
		cout << "sr_length: " << sr_length << "num_deletions: " << num_deletions << endl;
		return start_index + sr_length - num_deletions;
	}
	cout << "sr_length: " << sr_length << "num_deletions: " << num_deletions << endl;
	return start_index + sr_length;
}

int computeEndIndexAndFrequencies(uint32_t * const nuc, uint32_t * const cigar, vector<ncount_t> &freq, int start_index, int read_index)
{
	uint32_t n, c, position, starting_position, num_elements;
	int num_deletions = 0;
	int sr_length = nuc[read_index];
	uint32_t * sr_start = &nuc[1 + read_index];
	uint32_t * cigar_start = &cigar[1 + read_index];
	num_elements = (int)ceil(sr_length / NUC_PER_ELEM);
	position = (num_elements * NUC_PER_ELEM) - sr_length;
	starting_position = position;
	//cout << "vector size: " << freq.size() << endl;
	for (int i = starting_position; i < (int) num_elements * NUC_PER_ELEM; i++) {
		n = nucleotide(sr_start, i);
		c = nucleotide(cigar_start, i);
		//cout << "start index: " << start_index << " position " << position << endl;
		// if match or mismatch, increment the corresponding nucleotide count
		if (c == 0 || c == 3) {
			switch (n) {
			case 0:
				freq[start_index + position - starting_position].num_a++;
				break;
			case 1:
				freq[start_index + position - starting_position].num_c++;
				break;
			case 2:
				freq[start_index + position - starting_position].num_g++;
				break;
			case 3:
				freq[start_index + position - starting_position].num_t++;
				break;
			}
			position++;
		}
		// if insertion, don't increment the position and do nothing
		else if (c == 2) {
			/* do nothing */
		}
		// if deletion, just increment the position
		else if (c == 1) {
			position++;
			num_deletions++;
		}
		//cout << "A: " << freq[i].num_a << " C: " << freq[i].num_c << " G: " << freq[i].num_g << " T: " << freq[i].num_t << endl;
	}
		//cout << "sr_length: " << sr_length << "num_deletions: " << num_deletions << endl;
	return start_index + sr_length - num_deletions;
}

int checkRead(uint8_t * const threshold, int start_index, int end_index)
{
	cout << "start index: " << start_index << " end_index: " << end_index << endl;
	for (int i = start_index; i < end_index; i++) {
		if (threshold[i])
			return 1;
	}
	return 0;
}

int * createExtraArray(int k, int grain_size)
{
	int * extra = new int[3];
	extra[0] = (int)floor(k / 2);
	extra[2] = (int)ceil(k / 2);
	extra[1] = extra[0] + extra[2];
	return extra;
}