#include "mapping.h"

/* Creates a hashtable for read sequence. This uses the C++ Standard Library
		seq_len: total number of bytes bytes
*/
void createStdMap(uint32_t *read_seq, int k, multimap<uint64_t, h_pair> * hashmap, int seq_len) {
	if (k <= 0 || k > 16) {
		return;
	}
	int mask;
	h_pair hash_pair;
	uint64_t temp, map_elem;
	int leftindex, rightindex, leftmod, rightmod, leftbound, rightbound;
	int num_reads = seq_len / READ_SIZE;
	for (int j = 0; j < num_reads; j++) { //j is the read number
		for (int i = NUM_NUCL; i > (int)((NUM_NUCL) - (read_seq[j*READ_SIZE] - k + 1)); i--) { //i is the nucleotide number within a read

			rightbound = i + j * NUM_NUCL;
			leftbound = i - k + j * NUM_NUCL;
			leftindex = (leftbound-1) / 16;
			rightindex = (rightbound-1) / 16;
			leftmod = (leftbound - 1) % 16;
			rightmod = (rightbound - 1) % 16;
			map_elem = 0;
			for (int a = leftindex; a <= rightindex; a++) {
				if (a == leftindex && a == rightindex) {
					temp = read_seq[a] >> (32 - (rightmod + 1)* 2);
					mask = (int)(0x1 << (k * 2)) - 1;
					map_elem = temp & mask;
				}
				else if (a == leftindex) {
					if (!leftmod) {
						mask = 0xffffffff;
					}
					else {
						mask = (int)(0x1 << ((16 - leftmod) * 2)) - 1;
					}
					temp = read_seq[a] & mask;
					temp = temp << (rightmod * 2);
					map_elem = temp;
				}
				else if (a == rightindex) {
					mask = (int)(0x1 << (rightmod * 2)) - 1;
					temp = (read_seq[a] >> (rightmod * 2)) & mask; 
					map_elem += temp;
				}
				else {
					temp = read_seq[a];
				}
			}
			//cout << "map_elem " << hex << map_elem << endl;
			hash_pair = std::pair<uint32_t, uint32_t>(j, 128-i);
			hashmap->insert(std::pair<std::uint64_t, h_pair>(map_elem, hash_pair));
		}
	}
}

/* Creates a hashtable for read sequence. This does not use the C++ Standard Library
seq_len: total number of bytes bytes
*/
void createMap(uint32_t *read_seq, int k, uint32_t *table, int seq_len) {
	if (k <= 0 || k > 16) {
		return;
	}
	int mask, ind = 0;
	uint32_t temp, map_elem;
	int leftindex, rightindex, leftmod, rightmod, leftbound, rightbound;
	int num_reads = seq_len / READ_SIZE;
	for (int j = 0; j < num_reads; j++) { //j is the read number
		for (int i = NUM_NUCL; i >(int)((NUM_NUCL)-(read_seq[j*READ_SIZE] - k + 1)); i--) { //i is the nucleotide number within a read

			rightbound = i + j * NUM_NUCL;
			leftbound = i - k + j * NUM_NUCL;
			leftindex = (leftbound - 1) / 16;
			rightindex = (rightbound - 1) / 16;
			leftmod = (leftbound - 1) % 16;
			rightmod = (rightbound - 1) % 16;
			map_elem = 0;
			for (int a = leftindex; a <= rightindex; a++) {
				if (a == leftindex && a == rightindex) {
					temp = read_seq[a] >> (32 - (rightmod + 1) * 2);
					mask = (int)(0x1 << (k * 2)) - 1;
					map_elem = temp & mask;
				}
				else if (a == leftindex) {
					if (!leftmod) {
						mask = 0xffffffff;
					}
					else {
						mask = (int)(0x1 << ((16 - leftmod) * 2)) - 1;
					}
					temp = read_seq[a] & mask;
					temp = temp << (rightmod * 2);
					map_elem = temp;
				}
				else if (a == rightindex) {
					mask = (int)(0x1 << (rightmod * 2)) - 1;
					temp = (read_seq[a] >> (rightmod * 2)) & mask;
					map_elem += temp;
				}
				else {
					temp = read_seq[a];
				}
			}
			//cout << "map_elem " << hex << map_elem << endl;
			table[ind] = map_elem;
			ind++;
		}
	}
}

unsigned int nextPowerOf2(unsigned int n)
{
	unsigned int p = 1;
	if (n && !(n & (n - 1)))
		return n;

	while (p < n) {
		p <<= 1;
	}
	return p;
}

void countFrequency(uint32_t * table_in, kmer_t * table_out, int size)
{
	int out_ind = 0;
	table_out[out_ind] = kmer_t{table_in[out_ind], 1};
	for (int i = 1; i < size; i++) {
		if (table_out[out_ind].kmer == table_in[i]) {
			table_out[out_ind].k_num++;
		}
		else {
			out_ind++;
			table_out[out_ind].kmer = table_in[i];
			table_out[out_ind].k_num = 1;
		}
	}
}

/*int main(void)
{
	uint32_t * sample = createSampleReads(294912);
	//uint32_t * sample = createSmallSample();
	print_two_reads(sample);
	cout << endl;
	//multimap<uint64_t, h_pair> hashmap; 
	int k = 8;
	int num_kmers = (128 - k + 1)* 32 * 1024;
	//int num_kmers = (128 - k + 1) * 2;
	uint32_t *table = new uint32_t[num_kmers];
	auto start_time = std::chrono::high_resolution_clock::now();
	createMap(sample, k, table, LENGTH);
	auto current_time = std::chrono::high_resolution_clock::now();
	std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::duration<float>>(current_time - start_time).count() << " seconds" << std::endl;
	print_two_reads(table);
	system("pause");
	delete [] sample;
	delete[] table;
	return 0;
}*/