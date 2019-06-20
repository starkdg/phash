#include "pHash.h"

static int _bmb_setbit(BMBHash &bh, uint32_t bit) {
	if (bh.hash == NULL) return -1;
	uint8_t one = 0x01;
	bh.hash[bit/8] |= one << (bit%8);
	return 0;
}
static void _ph_bmb_new(BMBHash &bh, uint32_t bytelength) {
	bh.bytelength = bytelength;
	bh.hash = (uint8_t*)calloc(sizeof(uint8_t), bytelength);
}

void ph_bmb_free(BMBHash &bh) {
	free(bh.hash);
}

int ph_bmb_imagehash(const char *file, BMBHash &ret_hash) {
	if (!file) return -1;

	CImg<uint8_t> img;
	const int preset_size_x = 256;
	const int preset_size_y = 256;
	const int blk_size_x = 16;
	const int blk_size_y = 16;

	int number_of_blocks = preset_size_x/blk_size_x*preset_size_y/blk_size_y; 
	uint32_t bytelength = number_of_blocks/8;
	
	try {
		img.load(file);
	} catch (CImgIOException ex) {
		return -1;
	}

	switch (img.spectrum()) {
	case 3: // from RGB
		img = img.RGBtoYUV().channel(0);
		break;
	case 1: // grayscale
		break;;
	default:
		return -1;
	}
	img.resize(preset_size_x, preset_size_y);

	double *mean_vals = new double[number_of_blocks];

	int blockidx = 0;
	for (int blockrow = 0; blockrow < preset_size_y - blk_size_y; blockrow += blk_size_y) {
		for (int blockcol = 0; blockcol < preset_size_x - blk_size_x; blockcol += blk_size_x) {
			CImg<uint8_t> subimg = img.crop(blockcol, blockrow, 0, 0, blockcol+blk_size_x, blockrow+blk_size_y, 0, 0, 0);
			mean_vals[blockidx] = subimg.mean();
			blockidx++;
		}
	}

	/* calculate the median */
	double median_value = CImg<double>(mean_vals, number_of_blocks).median();
	
	_ph_bmb_new(ret_hash, bytelength);

	for (int i = 0; i < number_of_blocks; i++) {
		if (mean_vals[i] <= median_value) {
			_bmb_setbit(ret_hash, 0);
		} else {
			_bmb_setbit(ret_hash, 1);
		}
	}
	delete[] mean_vals;
   	return 0;
}

int ph_bmb_distance(const BMBHash &bh1, const BMBHash &bh2){

	return ph_hammingdistance2(bh1.hash, bh1.bytelength, bh2.hash, bh2.bytelength);
}

