
#include "parse_gen.h"
#include <unordered_map>
#include <stdexcept>
#include "gsl/gsl_cdf.h"
#include "math.h"
#include <fstream>
#include <sstream>

using std::ifstream; using std::string;
using std::vector; using std::unordered_map; 
using std::cout; using std::endl;
using std::pair;

double sign(double x) {
    if (x > 0) return 1.0;
    else if (x < 0) return -1.0;
    else return 0;
}

void parse_ref(const string &ref_path, \
	unordered_map<string, CoordInfo*> &ref_dict, \
	vector<pair<size_t, size_t>> &boundary, \
	vector<string> &SNP) {
    ifstream infile(ref_path.c_str());
    string id, A1, A2, line;

    if (!infile) {
	throw std::runtime_error("Error: cannot open ref "
		"snpInfo file: " + ref_path);
    }

    int n = 0, section = 1;
    size_t left, right;
   
    // skip the header
    getline(infile, line);
    while (getline(infile, line)) {
	if (line == "") {
	    section++;
	    cout << "Readed " << n << " LD blocks." << endl;
	    n = 0;
	    getline(infile, line); // skip the header
	    continue;
	}

	if (section == 1) {
	    std::istringstream my_stream(line);
	    my_stream >> left >> right;
	    boundary.push_back(std::make_pair(left, right));
	    n++;
	}
	else {
	    std::istringstream my_stream(line);
	    my_stream >> id >> A1 >> A2;
	    SNP.push_back(id);
	    CoordInfo *ref_info = new CoordInfo;
	    ref_info->A1 = A1; ref_info->A2 = A2;
	    ref_info->include_ref = false; 
	    ref_info->include_ss = false;
	    ref_info->beta = 0;
	    if (!ref_dict.insert(pair<string, \
			CoordInfo*>(id, ref_info)).second) {
		throw std::runtime_error("Error: duplicate SNP found "
			"in ref snpInfo file: " + id);
	    }
	    n++;
	}
    }

    cout << "Readed " << n << " SNPs from reference panel" << endl;

    infile.close();
}

void parse_valid(const string &valid_path, \
	unordered_map<string, CoordInfo*> &ref_dict) {
    ifstream infile(valid_path.c_str());
    if (!infile) {
	throw std::runtime_error("Error: cannot open "
		"ref snpInfo file: " + valid_path);
    }

    string id, A1, A2, header;
    float genPos; unsigned chr, phyPos;
    int n = 0;
    unordered_map<string, CoordInfo*>::iterator idx;

    while (infile >> chr >> id >> genPos >> phyPos >> A1 >> A2) {
	idx = ref_dict.find(id);
	if (idx != ref_dict.end()) {
	    idx->second->include_ref = true;
	    n++;
	}
    }
    cout << n << " common SNPs between reference "
	"and validation datasets." << endl;
    infile.close();
}

void parse_ss(const string &ss_path, unordered_map<string, \
	CoordInfo*> &ref_dict, unsigned sz) {
    ifstream infile(ss_path.c_str());
    if (!infile) {
	throw std::runtime_error("Error: cannot open "
		"summary statistics: " + ss_path);
    }

    string id, A1, A2, header;
    double beta, pval;
    int n = 0;
    getline(infile, header);
    unordered_map<string, CoordInfo*>::iterator idx;

    int n_flip = 0, n_bad = 0;
    while (infile >> id >> A1 >> A2 >> beta >> pval) {
	idx = ref_dict.find(id);
	if (pval <= 1e-323) {
	    pval = 1e-323;
	}
	if (idx != ref_dict.end() && idx->second->include_ref) {
	    if (A1 == idx->second->A1 && A2 == idx->second->A2) {
		idx->second->include_ss = true;
		idx->second->beta = 1.0*sign(beta)* \
		fabs(gsl_cdf_ugaussian_Pinv(pval/2.0))/sqrt(sz);
		n++;
	    }
	    else if (A1 == idx->second->A2 && A2 == idx->second->A1) {
		idx->second->include_ss = true;
		idx->second->beta = -1.0*sign(beta)* \
		fabs(gsl_cdf_ugaussian_Pinv(pval/2.0))/sqrt(sz);
		n++;
		n_flip++;
	    }
	    else {
		n_bad++;
	    }
	}
    }

    cout << n_flip << " SNPs have flipped alleles between summary statistics and " 
	<< "reference panel." << endl;
    cout << n_bad << " SNPs removed due to mismatch of allels between " 
	<< "summary statistics and reference panel." << endl;
    cout << n << " common SNPs among reference, validation "
	"and gwas summary statistics." << endl;
    infile.close();
}

void parse_ld_mat(const string &ldmat_path, unordered_map<string, CoordInfo*> &ref_dict, \
	const vector<pair<size_t, size_t>> &boundary, \
	const vector<string> &SNP, mcmc_data &dat) {
    
    FILE *fp;
    fp = fopen(ldmat_path.c_str(), "rb");

    if (!fp) {
	throw std::runtime_error("Error: cannot open LD matrix file: " + ldmat_path);
    }

    unordered_map<string, CoordInfo*>::iterator idx;
    vector<size_t> snp_idx;
    size_t left = 0, right = 0;
    for (size_t i=0; i<boundary.size(); i++) {
	snp_idx.clear();
	for (size_t j=boundary[i].first; j<boundary[i].second; j++) {
	    idx = ref_dict.find(SNP[j-boundary[0].first]);
	    if (idx->second->include_ss) {
		dat.id.push_back(SNP[j-boundary[0].first]);
		dat.A1.push_back(idx->second->A1);
		dat.A2.push_back(idx->second->A2);
		dat.beta_mrg.push_back(idx->second->beta);
		snp_idx.push_back(j-boundary[i].first);
		right++;
	    }
	}
	gsl_matrix *tmp_mat = gsl_matrix_alloc(boundary[i].second-boundary[i].first, \
		                boundary[i].second-boundary[i].first);
	gsl_matrix_fread(fp, tmp_mat);

	if (left == right) {
	    gsl_matrix_free(tmp_mat);
	    continue; 
	}
	dat.boundary.push_back(std::make_pair(left,right));
	gsl_matrix *tmp_mat_sub = gsl_matrix_alloc(right-left, right-left);

	// copy rows from original matrix to second matrix with correct SNPs 
	for (size_t j=0; j<snp_idx.size(); j++) {
	    for (size_t k=0; k<snp_idx.size(); k++) {
		double tmp = gsl_matrix_get(tmp_mat, snp_idx[j], snp_idx[k]);
		gsl_matrix_set(tmp_mat_sub, j, k, tmp);
	    }
	}
	dat.ref_ld_mat.push_back(tmp_mat_sub);
	left = right;
	gsl_matrix_free(tmp_mat);
    }

    fclose(fp);
}
    

void coord(const string &ref_path, const string &ss_path, \
	const string &valid_path, const string &ldmat_path, \
	mcmc_data &dat, unsigned sz) {
    unordered_map<string, CoordInfo*> ref_dict;
    vector<pair<size_t, size_t>> boundary;
    vector<string> SNP;
    unordered_map<string, CoordInfo*>::iterator it;

    parse_ref(ref_path, ref_dict, boundary, SNP);
    
    if (!valid_path.empty()) {
	parse_valid(valid_path, ref_dict);
    }
    else {
	for (it=ref_dict.begin(); it != ref_dict.end(); it++) {
	    it->second->include_ref = true;
	}
    }

    parse_ss(ss_path, ref_dict, sz);    
    parse_ld_mat(ldmat_path, \
	    ref_dict, boundary, SNP, dat);
    for (it=ref_dict.begin(); it != ref_dict.end(); it++) {
	 delete(it->second);
    }    
}

/*int main() {
    mcmc_data dat;
    coord("test.snpInfo", "/ysm-gpfs/pi/zhao-data//gz222/height/summ_stats/PRS_cs.txt", \
	    "/ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3.bim", \
	    "test.dat", dat, 252230);
    for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
	gsl_matrix_free(dat.ref_ld_mat[i]);
    }

    return 0;
}*/

