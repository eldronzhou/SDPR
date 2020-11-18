
#include "SDPR_io.h"
#include "math.h"
#include <algorithm>

using std::cout; using std::endl; using std::ifstream;
using std::vector; using std::string;
using std::map;

// does not check validatiy of fam
size_t get_nsamples(const string &fam_path) {
   size_t n_samples = 0;
   FILE *fam_fp;
   fam_fp = fopen(fam_path.c_str(), "r");
   
   if (fam_fp == NULL) {
       throw("Error: cannot open fam file: " + fam_path + ".");
   }

   char line[2048];
   while (fgets(line, 2048, fam_fp) != NULL)
       n_samples++;

   fclose(fam_fp);

   return n_samples;
}

// read genotypes into a vector from .bed file for n th marker (0-based)
// in .bim file
void read_bed(gsl_matrix *snp, const string &bed_path, \
	const size_t n_samples, size_t left, size_t right) {
    FILE *bed_fp;
    bed_fp = fopen(bed_path.c_str(), "rb");

    if (bed_fp == NULL) {
	throw("Error: cannot open bed file: " + bed_path + ".");
    }

    char header[3];
    fread(header, sizeof(header), 1, bed_fp);

    if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
	throw("Incorrect first three bytes of bed file: " + bed_path + ".");
    }
 
    size_t size;
    if (n_samples % 4 == 0) {
	size = n_samples / 4;
    }
    else {
	size = n_samples / 4 + 1;
    }

    fseek(bed_fp, left*size+3, SEEK_SET);

    size_t n_snps = right - left;

    double *geno = new double[n_samples];

    for (size_t snp_idx=0; snp_idx<n_snps; snp_idx++) {

	size_t n_miss = 0; double sum = 0.0;
	
	size_t idx = 0;

	// convert hex to binary, then to genotype
	// 00 (2) homo for first allele in .bim
	// 01 (-9) missing genotype
	// 10 (1) heterozygous
	// 11 (0) homozygous for second allele in .bim
	
	int count[3] = {0};
	for (size_t j=0; j<size; j++) {
	    int v = 0;
	    fread(&v, 1, 1, bed_fp);
	    
	    // genotype were read backwards
	    for (size_t i=0; i<4; i++) {
		if (4*j+i+1 > n_samples) continue;

		if (v % 2 == 0) {
		    v /= 2;
		    if (v % 2 == 0) {
			geno[idx] = 0;
			count[0]++;
		    }
		    else {
			geno[idx] = 1;
			sum++;
			count[1]++;
		    }
		}
		else {
		    v /= 2;
		    if (v % 2 == 0) {
			geno[idx] = -9;
			n_miss++;
		    }
		    else {
			geno[idx] = 2;
			sum += 2;
		    }
		}
		
		v /= 2;
		idx++;
	    }
	}

	// fill in missing genotype by the most frequent allele
	double mean = sum / (n_samples-n_miss);

	if (n_miss > 0) {
	    for (size_t i=0; i<n_samples; i++) {
		if (geno[i] == -9) {
		     //geno[i] = mean;
		    geno[i] = std::max_element(count, count+3)-count;
		}
	    }
	}

	// standardize the genotype
	mean = gsl_stats_mean(geno, 1, n_samples);
	
	double sd = gsl_stats_sd_m(geno, 1, n_samples, mean)*sqrt(n_samples-1)/sqrt(n_samples);

	if (sd == 0) sd = 1;

	for (size_t i=0; i<n_samples; i++) {
	    gsl_matrix_set(snp, snp_idx, i, (geno[i]-mean)/sd);
	}
    }

    delete[] geno;
    fclose(bed_fp);
}

// read .bim file
void read_bim(const string &bim_path, SnpInfo *snpinfo) {

    ifstream infile(bim_path.c_str());
    if (!infile) {
	throw("Error: cannot open bim file: " + bim_path + ".");
    }

    cout << "Reading bim file from: " + bim_path + "." << endl;
   
    string id, A1, A2;
    float genPos; unsigned chr, phyPos;
    while (infile >> chr >> id >> genPos >> phyPos >> A1 >> A2) {
	/*if (std::find(snpinfo->id.begin(), snpinfo->id.end(), id) != \
		snpinfo->id.end()) {
	    throw("Error: Duplicated SNP ID " + id + " found in bim file.");
	}*/
	
	snpinfo->id.push_back(id);
	snpinfo->A1.push_back(A1);
	snpinfo->A2.push_back(A2);
	snpinfo->Pos.push_back(phyPos);
	if (snpinfo->chr_idx.find(chr) != snpinfo->chr_idx.end()) {
	    snpinfo->chr_idx[chr]++;
	}
	else {
	    cout << "Only support chr 1-22, but found chr: " << chr << endl;
	}
    }

    for (size_t i=1; i<23; i++) {
	if (snpinfo->chr_idx[i] != 0) {
	    snpinfo->chr_idx[i] += snpinfo->chr_idx[i-1];
	}
    }

    infile.close();
}

