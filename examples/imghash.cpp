/*

    pHash, the open source perceptual hash library
    Copyright (C) 2009 Aetilius, Inc.
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Evan Klinger - eklinger@phash.org
    David Starkweather - dstarkweather@phash.org

*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "pHash.h"

using namespace std;

namespace fs=boost::filesystem;
namespace po=boost::program_options;

enum ImgHash {DCTHASH=1, RADIALHASH, MHASH, BMBHASH };

int read_filenames_from_dir(const fs::path &dname, vector<fs::path> &files);

double compare_dct_imagehash(const fs::path &file1, const fs::path &file2);

double compare_radial_imagehash(const fs::path &file1, const fs::path &file2);

double compare_mh_imagehash(const fs::path &file1, const fs::path &file2);

double compare_bmb_imagehash(const fs::path &file1, const fs::path &file2);

double compare_images(const fs::path &file1, const fs::path &file2, const ImgHash method);

int calc_inter_distances(const vector<fs::path> &files1, const vector<fs::path> &files2,
					vector<double> &distances, const ImgHash method);

int calc_intra_distances(const vector<fs::path> &files, vector<double> &distances, const ImgHash method);

int find_stats(const vector<double> &values, double &minval, double &maxval, double &mean, double &var);

int compute_histogram(const vector<double> &values, const double max_value, const int nbins, vector<double> &histogram);


typedef struct Args {
	string dir1;
	string dir2;
	string outfile;
	bool help;
	int method;
	int nbins;
} Args;

Args parse_options(int argc, char **argv){
	Args args;
	po::options_description descr("Available Options");
	descr.add_options()
		("help,h", "usage information")
		("dir1", po::value<string>(&args.dir1)->required(), "path to directory 1")
		("dir2", po::value<string>(&args.dir2)->required(), "path to directory 2")
		("output, o", po::value<string>(&args.outfile)->default_value("histogram.dat"),
		                   "filename to write histogram values")
		("method,m", po::value<int>(&args.method)->default_value(1),
		                         "Specify image hash method to use.\n"
                                  "Method:\n"
		                          "\t 1 - dct image hash\n" 
		                          "\t 2 - radial image hash\n"
                                  "\t 3 - mh image hash\n"
		                          "\t 4 - bmb image hash\n\n")
		("nbins,b", po::value<int>(&args.nbins)->default_value(64),
		                          "number of histogram bins");

	po::variables_map vm;
	
	try {
		po::store(po::parse_command_line(argc, argv, descr), vm);
		po::notify(vm); 
		
		if (vm.count("help")){
			cout << endl << descr << endl << endl;
			exit(0);
		}
		
	} catch (exception &ex){
		cout << "ERROR: " << ex.what() << endl;
		cout << endl << descr << endl << endl;
		exit(0);
	}

	return args;
}


int main(int argc, char **argv){
	Args args = parse_options(argc, argv);
	
	const ImgHash method = (ImgHash)args.method;
	const int nbins = args.nbins;
	const string outfile = args.outfile;
	
	const fs::path dirname1(args.dir1);
	const fs::path dirname2(args.dir2);

	cout << "compare files in two directories" << endl;
	cout << "dir: " << dirname1 << endl;
	cout << "dir: " << dirname2 << endl;
	
	vector<fs::path> files1;
	int n1 = read_filenames_from_dir(dirname1, files1);
	assert(n1 >= 0);
	
	cout << "dir: " << dirname1 << " (" << n1 << " files)" << endl;

	vector<fs::path> files2;
	int n2 = read_filenames_from_dir(dirname2, files2);
	assert(n2 >= 0);
	assert(n1 == n2);

	cout << "dir: " << dirname2 << " (" << n2 << " files)" << endl;
	cout << "Calculate distances between similar files ..." << endl;

	vector<double> interdistances;
	int rc = calc_inter_distances(files1, files2, interdistances, method);
	assert(rc == 0);

	cout << "Calculate distances between dissimilar files ..." << endl;

	vector<double> intradistances;
	rc = calc_intra_distances(files1, intradistances, method);
	assert(rc == 0);

	/* find maximum of distance range */
	double min_d, max_d, mean_d, var_d;
	find_stats(intradistances, min_d, max_d, mean_d, var_d);


	/* calculate histograms */
	vector<double> inter_histogram;
	compute_histogram(interdistances, max_d, nbins, inter_histogram);

	vector<double> intra_histogram;
	compute_histogram(intradistances, max_d, nbins, intra_histogram);


	/* write out histograms to output file */
	ofstream ostrm(outfile, ios::out|ios::trunc);
	for (int i=0;i<nbins;i++){
		ostrm << left << setw(10) << i << left << setw(10) << inter_histogram[i]
			  << left << setw(10) << intra_histogram[i] << endl;
	}
	
	cout << "Done" << endl;
    return 0;
}


/**
 * **************************************************************** 
 *
 *	             Auxiliary functions below this point.
 *
 *
 * ****************************************************************
 **/

int read_filenames_from_dir(const fs::path &dname, vector<fs::path> &files){

	int n_files = 0;
	try {
		if (fs::exists(dname) && fs::is_directory(dname)){

			for (auto &&entry : fs::directory_iterator(dname)){
				files.push_back(entry.path());
				n_files++;
			}
			sort(files.begin(), files.end());
		} else {
			n_files = -1;
		}
	} catch (fs::filesystem_error &ex){
		cout << "ERROR reading from " << dname << " : " << ex.what() << endl;;
		return -1;
	}
	
	return n_files;
}

double compare_dct_imagehash(const fs::path &file1, const fs::path &file2){

	ulong64 hash1;
	if (ph_dct_imagehash(file1.c_str(), hash1) < 0)
		return -1.0;

	ulong64 hash2;
	if (ph_dct_imagehash(file2.c_str(), hash2) < 0)
		return -1.0;
		
	int d = ph_hamming_distance(hash1, hash2);
	return (double)d;
}

double compare_radial_imagehash(const fs::path &file1, const fs::path &file2){
	const double sigma = 1.0, gamma = 1.0;
	const int n_angles = 180;
	Digest digest1, digest2;

	if (ph_image_digest(file1.c_str(), sigma, gamma, digest1, n_angles) < 0)
		return -1;

	if (ph_image_digest(file2.c_str(), sigma, gamma, digest2, n_angles) < 0)
		return -1;

	double pcc;
	if (ph_peakcrosscorr(digest1, digest2, pcc) < 0)
		return -1;

	return 1.0 - pcc;
}

double compare_mh_imagehash(const fs::path &file1, const fs::path &file2){
	const float alpha = 2.0, lvl = 1.0;
	int n1, n2;
	uint8_t *mh1 = ph_mh_imagehash(file1.c_str(), n1, alpha, lvl);
	if (mh1 == NULL)
		return -1;
	
	uint8_t *mh2 = ph_mh_imagehash(file2.c_str(), n2, alpha, lvl);
	if (mh2 == NULL)
		return -1;
	
	int d = ph_hammingdistance2(mh1, n1, mh2, n2);
	return (double)d;
}

double compare_bmb_imagehash(const fs::path &file1, const fs::path &file2){
	BMBHash bh1, bh2;

	if (ph_bmb_imagehash(file1.c_str(), bh1) < 0)
		return -1;
	
	if (ph_bmb_imagehash(file2.c_str(), bh2) < 0)
		return -1;

	int d = ph_bmb_distance(bh1, bh2);

	ph_bmb_free(bh1);
	ph_bmb_free(bh2);
	
	return (double)d;
}

double compare_images(const fs::path &file1, const fs::path &file2, const ImgHash method){
	double d;
	
	switch (method){
	case DCTHASH:
		d = compare_dct_imagehash(file1, file2);
		break;
	case RADIALHASH:
		d = compare_radial_imagehash(file1, file2);
		break;
	case MHASH:
		d = compare_mh_imagehash(file1, file2);
		break;
	case BMBHASH:
		d = compare_bmb_imagehash(file1, file2);
		break;
	default:
		d = -1;
	}
	return d;
}

int calc_inter_distances(const vector<fs::path> &files1, const vector<fs::path> &files2,
					vector<double> &distances, const ImgHash method){
	if (files1.size() != files2.size() || files1.size() <= 0)
		return -1;
	distances.clear();
	for (int i=0;i < files1.size();i++){
		double d = compare_images(files1[i], files2[i], method);
		if (d < 0)
			continue;
		distances.push_back(d);
	}
	
	return 0;
}

int calc_intra_distances(const vector<fs::path> &files, vector<double> &distances, const ImgHash method){

	distances.clear();
	int ilimit = (files.size() >= 10) ? 10 : files.size();
	for (int i=0;i < ilimit-1;i++){
		int jlimit = (files.size() >= i+10) ? i+10 : files.size();
		for (int j=i+1;j < jlimit;j++){
			double d = compare_images(files[i], files[j], method);
			if (d < 0)
				continue;
			distances.push_back(d);
		}
	}

	return 0;
}

int find_stats(const vector<double> &values, double &minval, double &maxval, double &mean, double &var){
	int n_values = 0;
	double sum = 0;
	double sumsq = 0;
	double min_v = numeric_limits<double>::max();
	double max_v = numeric_limits<double>::min();
	for (double val : values){
		sum = sum + val;
		sumsq = sumsq + val*val;
		n_values++;

		if (val > max_v)
			max_v = val;
		if (val < min_v)
			min_v = val;
	}

	mean /= n_values;
	var = sumsq - ((sum*sum)/n_values)/(n_values-1);
	minval = min_v;
	maxval = max_v;
	return 0;
}

int compute_histogram(const vector<double> &values, const double max_value, const int nbins, vector<double> &histogram){
	histogram.clear();
	histogram.resize(nbins, 0);

	/* count values by bin */
	double bin_width = max_value/(double)nbins;
	for (double val : values){
		int bin = (int)floor(val/bin_width);
		histogram[bin] += 1;
	}
	
	/* normalize histogram to number of values */
	for (int i=0;i < nbins;i++){
		histogram[i] /= (double)values.size();
	}
	
	return 0;
}

