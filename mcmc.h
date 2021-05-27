
#include <iostream>
#include <vector>
#include "parse_gen.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "math.h"

typedef struct {
    std::vector<gsl_matrix*> A;
    std::vector<gsl_matrix*> B;
    std::vector<gsl_matrix*> L;
    std::vector<gsl_vector*> beta_mrg;
    std::vector<gsl_vector*> calc_b_tmp;
    std::vector<double> num;
    std::vector<double> denom;
} ldmat_data;

class MCMC_state {
    public:
	double alpha;
	double eta;
	double h2;
	double N;
	gsl_vector *beta;
	gsl_vector *b;
	std::vector<int> cls_assgn;
	std::vector<double> V;
	std::vector<double> p;
	std::vector<double> log_p;
	std::vector<double> cluster_var;
	std::vector<unsigned> suff_stats;
	std::vector<double> sumsq;
	MCMC_state(size_t num_snp, size_t max_cluster, \
		double a0, double b0, double sz) {
	    a0k = a0; b0k = b0; N = sz;
	    // Changed May 20 2021
	    // Now N (sz) is absorbed into A, B; so set to 1.
	    N = 1.0;

	    n_snp = num_snp;
	    M = max_cluster;
	    alpha = 1;
	    eta = 1;
	    beta = gsl_vector_calloc(num_snp);
	    b = gsl_vector_calloc(num_snp);
	    p.assign(max_cluster, 1.0/max_cluster);
	    log_p.assign(max_cluster, 0);
	    for (size_t i=0; i<max_cluster; i++) {
		log_p[i] = logf(p[i] + 1e-40);
	    }

	    cluster_var.assign(max_cluster, 0.0);
	    suff_stats.assign(max_cluster, 0);
	    sumsq.assign(max_cluster, 0.0);
	    V.assign(max_cluster, 0.0);
	    cls_assgn.assign(num_snp, 0);
	    r = gsl_rng_alloc(gsl_rng_default);
	    for (size_t i=0; i<num_snp; i++) {
		cls_assgn[i] = gsl_rng_uniform_int(r, M);
	    }
	}

	~MCMC_state() {
	    gsl_vector_free(beta);
	    gsl_vector_free(b);
	    gsl_rng_free(r);
	}

	void sample_sigma2();
	void calc_b(size_t j, const mcmc_data &dat, const ldmat_data &ldmat_dat);
	void sample_assignment(size_t j, const mcmc_data &dat, \
		        const ldmat_data &ldmat_dat);
	void update_suffstats();
	void sample_V();
	void update_p();
	void sample_alpha();
	void sample_beta(size_t j, const mcmc_data &dat, \
		       ldmat_data &ldmat_dat);
	void compute_h2(const mcmc_data &dat);
	void sample_eta(const ldmat_data &ldmat_dat);

    private:
	double a0k;
	double b0k;
	size_t M, n_snp;
	gsl_rng *r;
};

class MCMC_samples {
    public:
	gsl_vector *beta;
	double h2;

	MCMC_samples(size_t num_snps) {
	    beta = gsl_vector_calloc(num_snps);
	    h2 = 0;
	}

	~MCMC_samples() {
	    gsl_vector_free(beta);
	}
};

void mcmc(const std::string &ref_path, const std::string &ss_path, \
	const std::string &valid_path, const std::string &ldmat_path, \
	const std::string &out_path, unsigned sz, double a, double c, \
	size_t M, double a0k, double b0k, \
	int iter, int burn, int thin, unsigned n_threads, int opt_llk);



