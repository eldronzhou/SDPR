
#include "SDPR_io.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include <thread>
#include <algorithm>

using std::cout; using std::endl;
using std::string; using std::vector;
using std::pair; using std::ofstream; 
using std::thread; using std::ref;

void calc_ref_ld_shrink(size_t k, gsl_matrix **ref_ld_mat, const string &bed_path, \
	const vector <pair<size_t, size_t>> &boundary, size_t n_sample) {
    size_t left = boundary[k].first;
    size_t right = boundary[k].second;
    gsl_matrix *snp = gsl_matrix_calloc(right-left, n_sample);

    read_bed(snp, bed_path, n_sample, left, right);

    gsl_matrix *snp2 = gsl_matrix_calloc(right-left, n_sample);
    gsl_matrix_memcpy(snp2, snp);
    gsl_matrix_mul_elements(snp2, snp2);

    gsl_matrix *snp2_prod = gsl_matrix_calloc(right-left, right-left);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, snp2, snp2, 0.0, snp2_prod);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, snp, snp, 0.0, ref_ld_mat[k]);
    gsl_matrix_scale(ref_ld_mat[k], 1.0/n_sample);

    double num = 0, denom = 0;

    for (size_t i = 0; i<right-left; i++) {
	for (size_t j=0; j<right-left; j++) {
	    if (i == j) continue;
	    num +=  gsl_matrix_get(snp2_prod, i, j) * n_sample / gsl_pow_3((double) n_sample-1) - \
	     gsl_pow_2(gsl_matrix_get(ref_ld_mat[k], i, j)) * gsl_pow_2((double) n_sample) / gsl_pow_3((double) n_sample-1); 
	    
	    denom += gsl_pow_2(gsl_matrix_get(ref_ld_mat[k], i, j));
	}
    }

    double sr = 1.0 - num/denom;
    if (sr < 0) sr = 0;
    if (sr > 1) sr = 1;

    gsl_matrix *tg_mat = gsl_matrix_alloc(right-left, right-left);
    gsl_matrix_set_identity(tg_mat);
    gsl_matrix_scale(ref_ld_mat[k], sr);
    gsl_matrix_scale(tg_mat, 1.0-sr);

    gsl_matrix_add(ref_ld_mat[k], tg_mat);

    gsl_matrix_free(snp); gsl_matrix_free(snp2);
    gsl_matrix_free(snp2_prod); gsl_matrix_free(tg_mat);
}

void calc_ref_parallel(size_t i, const vector<size_t> *v, gsl_matrix **ref_ld_mat, \
	const string &bed_path, const vector <pair<size_t, size_t>> &boundary, size_t n_sample) {

    for (size_t k=0; k<v[i].size(); k++) {
	calc_ref_ld_shrink(v[i][k], ref_ld_mat, bed_path, boundary, n_sample);
    }
}

bool myCmp(const pair<size_t, size_t> &a, const pair<size_t, size_t> &b) {
    return a.second > b.second;
}

void div_block(const string &pfile, \
	const string &out_dir, \
	unsigned chrom, size_t n_thread, double r2) {
    string fam_path = pfile + ".fam";
    string bim_path = pfile + ".bim";
    string bed_path = pfile + ".bed";
    size_t n_sample = get_nsamples(fam_path.c_str());
    SnpInfo snpinfo;
    for (size_t i=0; i<23; i++) {
	snpinfo.chr_idx[i] = 0;
    }
    read_bim(bim_path.c_str(), &snpinfo);

    for (size_t i=0; i<23; i++) {
	cout << "chrom " << i+1 << " " << snpinfo.chr_idx[i]  << endl;
    }

    size_t left = snpinfo.chr_idx[chrom-1], right = snpinfo.chr_idx[chrom];

    gsl_matrix *snp = gsl_matrix_calloc(right-left, n_sample);

    read_bed(snp, bed_path, n_sample, left, right);

    cout << "Readed " << right - left << " SNPs on Chr " << chrom << endl;

    // divide into approx. indep. blocks
    double cor = 0; 
    size_t *max_list = new size_t[right-left];
    gsl_vector_view snp1, snp2;
    for (size_t i=left; i<right; i++) {
	max_list[i-left] = i;
	snp1 = gsl_matrix_row(snp, i-left);

	for (size_t j=i+1; j<i+300; j++) {
	    if (j >= right) continue;
	    snp2 = gsl_matrix_row(snp, j-left);
	    gsl_blas_ddot(&snp1.vector, &snp2.vector, &cor);
	    cor /= n_sample;

	    if (cor*cor > r2) {
		max_list[i-left]= j;
	    }
    	}

	if (i == left) continue;

	if (max_list[i-left] < max_list[i-left-1]) {
	    max_list[i-left] = max_list[i-left-1];
	}
    }

    gsl_matrix_free(snp);

    vector<pair<size_t, size_t>> boundary;
    vector<pair<size_t, size_t>> blk_size;

    size_t left_bound = left, n_blk = 0;
    for (size_t i=left; i<right; i++) {
	if (max_list[i-left] == i) {
	    if (i+1-left_bound < 300 && i != right-1) continue;
	    boundary.push_back(std::make_pair(left_bound,i+1));
	    blk_size.push_back(std::make_pair(n_blk, i+1-left_bound));
	    left_bound = i+1;
	    n_blk++; 
	}
    }

    std::sort(blk_size.begin(), blk_size.end(), myCmp);

    cout << "Divided into " << n_blk << " indpenent blocks with max size: " \
	<< blk_size[0].second << endl;

    // calculate shrinkage ref ld mat

    gsl_matrix **ref_ld_mat = new gsl_matrix*[n_blk];
    for (size_t i=0; i<n_blk; i++) {
	ref_ld_mat[i] = gsl_matrix_calloc(boundary[i].second-boundary[i].first,\
		boundary[i].second-boundary[i].first);
    }

    vector<thread> threads(n_thread);
    
    unsigned *bin = new unsigned[n_thread];

    for (size_t i=0; i<n_thread; i++) {
	bin[i] = 0;
    }

    vector<size_t> *v = new vector<size_t>[n_thread];
    
    // binpacking to assign workitems to threads
    for (size_t i=0; i<n_blk; i++) {
	size_t idx = std::min_element(bin, bin+n_thread) - bin;
	bin[idx] += blk_size[i].second*blk_size[i].second;
	v[idx].push_back(blk_size[i].first);
    }

    /*for (size_t j=0; j<n_thread; j++) {
	cout << bin[j] << " ";
    }*/

    for (size_t i=0; i<n_thread; i++) {
	threads[i] = thread(calc_ref_parallel, i, ref(v), ref(ref_ld_mat), bed_path, ref(boundary), n_sample);
    }
    
    for (size_t i=0; i<n_thread; i++) {
	threads[i].join();
    }

    string out_ldmat = out_dir + "/chr" + \
		       std::to_string(chrom) + ".dat";
    FILE *f = fopen(out_ldmat.c_str(), "wb");
    for (size_t i=0; i<n_blk; i++) {
	gsl_matrix_fwrite(f, ref_ld_mat[i]);
	gsl_matrix_free(ref_ld_mat[i]);
    }
    fclose(f);
    
    string out_snpinfo = out_dir + "/chr" + \
			 std::to_string(chrom) + ".snpInfo";
    ofstream out(out_snpinfo);

    out << "start" << '\t' << "end" << endl;

    for (size_t i=0; i<boundary.size(); i++) {
	out << boundary[i].first << '\t' << boundary[i].second << endl;
    }

    out << endl;

    out << "SNP" << '\t' << "A1" << '\t' << "A2" << endl;
    for (size_t i=left; i<right; i++) {
	out << snpinfo.id[i] << '\t' << snpinfo.A1[i] \
	    << '\t' << snpinfo.A2[i] << endl;
    }
    out.close();

    delete[] ref_ld_mat; delete[] max_list; 
    delete[] bin; 
    delete[] v;
}

/*int main(int argc, char *argv[]) {

    div_block("tmp", 1, 4);

    return 0;
>*/


