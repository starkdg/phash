#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv){
	assert(argc == 3);
	const char* image_file1 = argv[1];
	const char* image_file2 = argv[2];

	assert(image_file1 != NULL);
	assert(image_file2 != NULL);

	ulong64 hash1;
	assert(ph_dct_imagehash(image_file1, hash1) == 1);

	ulong64 hash2;
	assert(ph_dct_imagehash(image_file2, hash2) == 1);
	
	int d = ph_hamming_distance(hash1, hash2);
	assert(d == 0);

	return 0;
}
