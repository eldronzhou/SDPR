
#include <algorithm>
#include "mcmc.h"
#include "math.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"
#include <thread>
#include <chrono>
#include "sse_mathfun.h"
#include "function_pool.h"
#include <fstream>
#include <x86intrin.h>
#include <numeric>

using namespace std::chrono;

using std::cout; using std::endl;
using std::thread; using std::ref;
using std::vector; using std::ofstream;
using std::string; using std::min;

#define square(x) ((x)*(x))

void MCMC_state::sample_sigma2() {
    for (size_t i=1; i<M; i++) {
	double a = suff_stats[i] / 2.0 + a0k;
	double b = 1.0 / (sumsq[i] / 2.0 + b0k);
	cluster_var[i] = 1.0/gsl_ran_gamma(r, a, b);
	if (isinf(cluster_var[i])) {
	    cluster_var[i] = 1e5;
	    std::cerr << "Cluster variance is infintie." << std::endl;
	}
	else if (cluster_var[i] == 0) {
	    cluster_var[i] = 1e-10;
	    std::cerr << "Cluster variance is zero." << std::endl;
	}
    }
}
	
void MCMC_state::calc_b(size_t j, const mcmc_data &dat, const ldmat_data &ldmat_dat) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i = dat.boundary[j].second;
    gsl_vector_view b_j = gsl_vector_subvector(b, start_i, end_i-start_i);
    gsl_vector_view beta_j = gsl_vector_subvector(beta, \
	    start_i, end_i-start_i);

    gsl_vector_const_view diag = gsl_matrix_const_diagonal(ldmat_dat.B[j]);
    
    // diag(B) \times beta
    gsl_vector_memcpy(&b_j.vector, &beta_j.vector);
    gsl_vector_mul(&b_j.vector, &diag.vector);

    // eta^2 * (diag(B) \times beta) - eta^2 * B beta
    gsl_blas_dsymv(CblasUpper, -eta*eta, ldmat_dat.B[j], &beta_j.vector, \
	    eta*eta, &b_j.vector);

    // eta^2 * (diag(B) \times beta) - eta^2 * B beta + eta * A^T beta_mrg
    gsl_blas_daxpy(eta, ldmat_dat.calc_b_tmp[j], &b_j.vector);
}

void MCMC_state::sample_assignment(size_t j, const mcmc_data &dat, \
	const ldmat_data &ldmat_dat) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i = dat.boundary[j].second;
   
    float **prob = new float*[end_i-start_i];
    float **tmp = new float*[end_i-start_i];
    vector<float> Bjj(end_i-start_i);
    vector<float> bj(end_i-start_i); 
    vector<float> rnd(end_i-start_i);

    float max_elem, log_exp_sum = 0;
    
    for (size_t i=0; i<end_i-start_i; i++) {
	prob[i] = new float[M];
	tmp[i] = new float[M];

	Bjj[i] = gsl_matrix_get(ldmat_dat.B[j], i, i);
	bj[i] = gsl_vector_get(b, start_i+i);
	
	prob[i][0] = log_p[0];
	rnd[i] = gsl_rng_uniform(r);
    }

    // N = 1.0 after May 21 2021
    float C = pow(eta, 2.0) * N;

    // auto vectorized
    for (size_t i=0; i<end_i-start_i; i++) {
	for (size_t k=1; k<M; k++) {
	    //cluster_var[k] = k*.5;
	    prob[i][k] = C * Bjj[i] * cluster_var[k] + 1;
	}
    }

    // unable to auto vectorize due to log
    // explicitly using SSE 
    __m128 _v, _m;
    for (size_t i=0; i<end_i-start_i; i++) {
	size_t k = 1;
	for (; k<M; k+=4) { // require M >= 4
	    _v = log_ps(_mm_loadu_ps(&prob[i][k]));
	    _mm_storeu_ps(&tmp[i][k], _v);
	}

	for (; k<M; k++) {
	    tmp[i][k] = logf(prob[i][k]);
	}
    }

    // auto vectorized
    for (size_t i=0; i<end_i-start_i; i++) {
	for (size_t k=1; k<M; k++) {
	    prob[i][k] = -0.5*tmp[i][k] + log_p[k] + \
			 square(N*bj[i]) * cluster_var[k] / (2*prob[i][k]);
	}
    }
    
    // non-vectorized version
    /*for (size_t i=0; i<end_i-start_i; i++) {
	for (size_t k=1; k<M; k++) {
	    prob[i][k] = -0.5*log( pow(eta, 2.0) * N * Bjj[i] * cluster_var[k] + 1) \
		    + log(p[k] + 1e-40) + pow(N * bj[i], 2.0) / \
		    (2 * (N * Bjj[i] * pow(eta, 2.0) + 1.0 / cluster_var[k]));
	}
    }*/


    for (size_t i=0; i<end_i-start_i; i++) {
	// non SSE version to find max
	//max_elem = *std::max_element(&prob[i][0], &prob[i][M-1]);

	// SSE version to find max
	// https://shybovycha.github.io/2017/02/21/speeding-up-algorithms-with-sse.html
	_v = _mm_loadu_ps(prob[i]);
	size_t k = 4;
	for (; k<M; k+=4) {
	    _v = _mm_max_ps(_v, _mm_loadu_ps(&prob[i][k]));
	}
	
	for (size_t m=0; m<3; m++) {
	    _v = _mm_max_ps(_v, _mm_shuffle_ps(_v, _v, 0x93));
	}

	_mm_store_ss(&max_elem, _v);

	for (; k<M; k++) {
	    max_elem = (max_elem > prob[i][k]) ? \
		    (max_elem) : (prob[i][k]);
	}

	// non SSE version log exp sum
	/*for (size_t k=0; k<M; k++) {
	    log_exp_sum += expf(prob[i][k] - max_elem);
	}*/

	// SSE version log exp sum
	_m = _mm_load1_ps(&max_elem);
	_v = exp_ps(_mm_sub_ps(_mm_loadu_ps(prob[i]), _m));
	
	k = 4;
	for (; k<M; k+=4) {
	    _v = _mm_add_ps(_v, \
		    exp_ps(_mm_sub_ps(_mm_loadu_ps(&prob[i][k]), _m)));
	}

	_v = _mm_hadd_ps(_v, _v);
	_v = _mm_hadd_ps(_v, _v);
	_mm_store_ss(&log_exp_sum, _v);

	for (; k<M; k++) {
	    log_exp_sum += expf(prob[i][k] - max_elem);
	}

	log_exp_sum = max_elem + logf(log_exp_sum);

	/*for (size_t k=0; k<M; k++) {
	    if (i == 0 && j == 0) cout << prob[i][k] << " ";
	}*/

	cls_assgn[i+start_i] = M-1;
	for (size_t k=0; k<M-1; k++) {
	    rnd[i] -= expf(prob[i][k] - log_exp_sum);
	    if (rnd[i] < 0) {
		cls_assgn[i+start_i] = k;
		break;
	    }
	}

	delete[] prob[i]; delete tmp[i];
    }
    delete[] prob; delete[] tmp;
}

void MCMC_state::update_suffstats() {
    std::fill(suff_stats.begin(), suff_stats.end(), 0.0);
    std::fill(sumsq.begin(), sumsq.end(), 0.0);
    for (size_t i=0; i<n_snp; i++) {
	suff_stats[cls_assgn[i]]++;
	double tmp = gsl_vector_get(beta, i);
	sumsq[cls_assgn[i]] += square(tmp);
    }
}

void MCMC_state::sample_V() {
    vector<double> a(M-1);

    a[M-2] = suff_stats[M-1];
    for (int i=M-3; i>=0; i--) {
	a[i] = suff_stats[i+1] + a[i+1];
    }

    for (size_t i=0; i<M-1; i++) {
	V[i] = gsl_ran_beta(r, \
		1 + suff_stats[i], \
		alpha + a[i]);
    }
    V[M-1] = 1;
}

void MCMC_state::update_p() {
    vector<double> cumprod(M-1);
    
    cumprod[0] = 1 - V[0];

    for (size_t i=1; i<M-1; i++) {
	cumprod[i] = cumprod[i-1] * (1 - V[i]);
	
	if (V[i] == 1) {
	    std::fill(cumprod.begin()+i+1, cumprod.end(), 0.0);
	    break;
	}
    }

    p[0] = V[0]; 
    for (size_t i=1; i<M-1; i++) {
	p[i] = cumprod[i-1] * V[i];
    }

    double sum = std::accumulate(p.begin(), p.end()-1, 0.0);
    if (1 - sum > 0) {
	p[M-1] = 1 - sum;
    }
    else {
	p[M-1] = 0;
    }

    for (size_t i=0; i<M; i++) {
	log_p[i] = logf(p[i] + 1e-40); 
    }
}

void MCMC_state::sample_alpha() {
    double sum = 0, m = 0;
    for (size_t i=0; i<M; i++) {
	if (V[i] != 1) {
	    sum += log(1 - V[i]);
	    m++;
	}
    }
    
    if (m == 0) m = 1;

    alpha = gsl_ran_gamma(r, .1+m-1, 1.0/(.1-sum));
}

void MCMC_state::sample_beta(size_t j, const mcmc_data &dat, \
	        ldmat_data &ldmat_dat) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i = dat.boundary[j].second;

    vector <size_t>causal_list;
    for (size_t i=start_i; i<end_i; i++) {
	if (cls_assgn[i] != 0) {
	    causal_list.push_back(i);
	}
    }

    gsl_vector_view beta_j = gsl_vector_subvector(beta, \
	                            start_i, end_i-start_i);

    gsl_vector_set_zero(&beta_j.vector);

    if (causal_list.size() == 0) {
	ldmat_dat.num[j] = 0;
	ldmat_dat.denom[j] = 0;
	return;
    }
    else if (causal_list.size() == 1) {
	double var_k = cluster_var[cls_assgn[causal_list[0]]];
	double bj = gsl_vector_get(b, causal_list[0]);
	double Bjj = gsl_matrix_get(ldmat_dat.B[j], \
		causal_list[0]-start_i, \
		causal_list[0]-start_i);
	double C = var_k / (N*var_k*square(eta)*Bjj + 1.0);
	double rv = sqrt(C)*gsl_ran_ugaussian(r) + C*N*bj;
	gsl_vector_set(&beta_j.vector, causal_list[0]-start_i, \
		rv);
	ldmat_dat.num[j] = bj*rv;
	ldmat_dat.denom[j] = square(rv)*Bjj;
	return;
    }

    gsl_vector *A_vec = gsl_vector_alloc(causal_list.size());

    gsl_vector *A_vec2 = gsl_vector_alloc(causal_list.size());

    double C = square(eta)*N; 

    double *ptr = (double*) malloc(causal_list.size() * \
	    causal_list.size() * sizeof(ptr));

    if (!ptr) {
	std::cerr << "Malloc failed for block " 
	    " may due to not enough memory." << endl;
    }

    for (size_t i=0; i<causal_list.size(); i++) {
	// N = 1.0 after May 21 2021
	// A_vec = N A[,idx].T beta_mrg = N A beta_mrg[idx]
	gsl_vector_set(A_vec, i, N*eta*gsl_vector_get(ldmat_dat.calc_b_tmp[j], \
		    causal_list[i]-start_i));
	
	// B = N B_gamma + \Sigma_0^{-1}
	// auto-vectorized
	for (size_t k=0; k<causal_list.size(); k++) {
	    if (i != k) {
		ptr[i*causal_list.size()+k] = C * \
		ldmat_dat.B[j]->data[ldmat_dat.B[j]->tda * \
		(causal_list[i]-start_i) + \
		causal_list[k]-start_i];
	    }
	    else {
		ptr[i*causal_list.size()+k] = C * \
		ldmat_dat.B[j]->data[ldmat_dat.B[j]->tda * \
		(causal_list[i]-start_i) + \
		causal_list[i]-start_i] + \
		1.0/cluster_var[cls_assgn[causal_list[i]]];
	    }
	   // gsl_matrix_set(B, i, k, tmp);
	}
    }

    gsl_vector_memcpy(A_vec2, A_vec);

    gsl_matrix_view B = gsl_matrix_view_array(ptr, \
	    causal_list.size(), causal_list.size());

    gsl_vector *beta_c = gsl_vector_alloc(causal_list.size());
    
    for (size_t i=0; i<causal_list.size(); i++) {
	gsl_vector_set(beta_c, i, gsl_ran_ugaussian(r));
    }

    // (N B_gamma + \Sigma_0^-1) = L L^T
    gsl_linalg_cholesky_decomp1(&B.matrix);

    // \mu = L^{-1} A_vec
    gsl_blas_dtrsv(CblasLower, CblasNoTrans, \
	    CblasNonUnit, &B.matrix, A_vec);

    // N(\mu, I)
    gsl_blas_daxpy(1.0, A_vec, beta_c);

    // X ~ N(\mu, I), L^{-T} X ~ N( L^{-T} \mu, (L L^T)^{-1} )
    gsl_blas_dtrsv(CblasLower, CblasTrans, \
	    CblasNonUnit, &B.matrix, beta_c);

    // compute eta related terms
    for (size_t i=0; i<causal_list.size(); i++) {
	gsl_matrix_set(&B.matrix, i, i, \
	C*gsl_matrix_get(ldmat_dat.B[j], 
	    causal_list[i]-start_i, \
	    causal_list[i]-start_i));
    }

    gsl_blas_ddot(A_vec2, beta_c, &ldmat_dat.num[j]);
    gsl_blas_dsymv(CblasUpper, 1.0, &B.matrix, \
	    beta_c, 0, A_vec);
    gsl_blas_ddot(beta_c, A_vec, &ldmat_dat.denom[j]);
    ldmat_dat.denom[j] /= square(eta);
    ldmat_dat.num[j] /= eta;

    for (size_t i=0; i<causal_list.size(); i++) {
	gsl_vector_set(&beta_j.vector, causal_list[i]-start_i, \
		gsl_vector_get(beta_c, i));
    } 

    gsl_vector_free(A_vec);
    gsl_vector_free(A_vec2);
    gsl_vector_free(beta_c);
    free(ptr);
}

void MCMC_state::compute_h2(const mcmc_data &dat) {

    double h2_tmp = 0;
    h2 = 0;
    for (size_t j=0; j<dat.ref_ld_mat.size(); j++) {
	size_t start_i = dat.boundary[j].first;
	size_t end_i = dat.boundary[j].second;
	gsl_vector *tmp = gsl_vector_alloc(end_i-start_i);
	gsl_vector_view beta_j = gsl_vector_subvector(beta, \
	    start_i, end_i-start_i);
	gsl_blas_dsymv(CblasUpper, 1.0, \
	    dat.ref_ld_mat[j], &beta_j.vector, 0, tmp);
	gsl_blas_ddot(tmp, &beta_j.vector, &h2_tmp);
	h2 += h2_tmp;
    }
}
									
void MCMC_state::sample_eta(const ldmat_data &ldmat_dat) {
    double num_sum = std::accumulate(ldmat_dat.num.begin(), \
	    ldmat_dat.num.end(), 0.0);

    double denom_sum = std::accumulate(ldmat_dat.denom.begin(), \
	    ldmat_dat.denom.end(), 0.0);

    denom_sum += 1e-6;

    eta = gsl_ran_ugaussian(r) * sqrt(1.0/denom_sum) + \
	  num_sum / denom_sum;
}

void solve_ldmat(const mcmc_data &dat, ldmat_data &ldmat_dat, \
	const double a, unsigned sz, int opt_llk) {
    for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
	size_t size = dat.boundary[i].second - dat.boundary[i].first;
	gsl_matrix *A = gsl_matrix_alloc(size, size);
	gsl_matrix *B = gsl_matrix_alloc(size, size);
	gsl_matrix *L = gsl_matrix_alloc(size, size);
	gsl_matrix_memcpy(A, dat.ref_ld_mat[i]);
	gsl_matrix_memcpy(B, dat.ref_ld_mat[i]);
	gsl_matrix_memcpy(L, dat.ref_ld_mat[i]);


	if (opt_llk == 1) {
	    // (R + aNI) / N A = R via cholesky decomp
	    // Changed May 21 2021 to divide by N
	    // replace aN with a
	    gsl_vector_view diag = gsl_matrix_diagonal(B);
	    gsl_vector_add_constant(&diag.vector, a);
	}
	else {
	    // R_ij N_s,ij / N_i N_j
	    // Added May 24 2021
	    for (size_t j=0; j<size ; j++) {
		for (size_t k=0; k<size; k++) {
		    double tmp = gsl_matrix_get(B, j, k);
		    // if genotyped on two different arrays, N_s = 0
		    size_t idx1 = j + dat.boundary[i].first;
		    size_t idx2 = k + dat.boundary[i].first;
		    if ( (dat.array[idx1] == 1 && dat.array[idx2] == 2) || \
			    (dat.array[idx1] == 2 && dat.array[idx2] == 1) ) {
			tmp = 0;
		    }
		    else {
			tmp *= min(dat.sz[idx1], dat.sz[idx2]) / \
			       (1.1 * dat.sz[idx1] * dat.sz[idx2]);
		    }
		    gsl_matrix_set(B, j, k, tmp);
		}
	    }

	    // force positive definite
	    // B = Q \Lambda Q^T
	    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(size);
	    gsl_matrix *evac = gsl_matrix_alloc(size, size);
	    gsl_matrix *eval = gsl_matrix_calloc(size, size);
	    gsl_vector_view eval_diag = gsl_matrix_diagonal(eval);
	    gsl_eigen_symmv(B, &eval_diag.vector, evac, w);

	    // get minium of eigen value
	    double eval_min = gsl_matrix_get(eval, 0, 0);
	    for (size_t k=1; k<size; k++) {
		double eval_k = gsl_matrix_get(eval, k, k);
		if (eval_k <= eval_min) {
		    eval_min = eval_k;
		}
	    }

	    // restore lower half of B
	    for (size_t j=0; j<size; j++) {
		for (size_t k=0; k<j; k++) {
		    double tmp = gsl_matrix_get(B, k, j);
		    gsl_matrix_set(B, j ,k, tmp);
		}
	    }

	    // if min eigen value < 0, add -1.1 * eval to diagonal
	    for (size_t j=0; j<size; j++) {
		if (eval_min < 0) {
		    gsl_matrix_set(B, j, j, \
			    1.0/dat.sz[j+dat.boundary[i].first] - 1.1*eval_min);
		}
		else {
		    gsl_matrix_set(B, j, j, \
			    1.0/dat.sz[j+dat.boundary[i].first]);
		}
	    }

	    gsl_matrix_free(evac);
	    gsl_matrix_free(eval);
	    gsl_eigen_symmv_free(w);
	}

	gsl_linalg_cholesky_decomp1(B);
	gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, \
		CblasNonUnit, 1.0, B, A);
	gsl_blas_dtrsm(CblasLeft, CblasLower, CblasTrans, \
		                CblasNonUnit, 1.0, B, A);

	// Changed May 21 2021 to divide by N 
	if (opt_llk == 1) {
	    gsl_matrix_scale(A, sz);
	}

	// B = RA
	// Changed May 21 2021 as A may not be symmetric
	//gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, L, A, 0, B);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, L, A, 0, B);
	
	// L = R %*% R;
	gsl_matrix_mul_elements(L, L);

	// memory allocation for A^T beta_mrg
	// Changed May 21 2021 from A to A^T
	gsl_vector *beta_mrg = gsl_vector_alloc(size);
	for (size_t j=0; j<size; j++) {
	    gsl_vector_set(beta_mrg, j, dat.beta_mrg[j+dat.boundary[i].first]);
	}
	gsl_vector *b_tmp = gsl_vector_alloc(size);

	//gsl_blas_dsymv(CblasUpper, 1.0, A, beta_mrg, 0, b_tmp);
	// Changed May 21 2021 from A to A^T
	gsl_blas_dgemv(CblasTrans, 1.0, A, beta_mrg, 0, b_tmp);

	ldmat_dat.A.push_back(A);
	ldmat_dat.B.push_back(B);
	ldmat_dat.L.push_back(L);
	ldmat_dat.calc_b_tmp.push_back(b_tmp);
	ldmat_dat.beta_mrg.push_back(beta_mrg);
	ldmat_dat.denom.push_back(0);
	ldmat_dat.num.push_back(0);
    }
}


void mcmc(const string &ref_path, const string &ss_path, \
	const string &valid_path, const string &ldmat_path, \
	const string &out_path, unsigned sz, double a, double c, \
	size_t M, double a0k, double b0k, \
	int iter, int burn, int thin, unsigned n_threads, int opt_llk) {

    // number of pst. samples 
    int n_pst = (iter-burn) / thin;

    cout << "Running SDPR with opt_llk " << opt_llk  << endl;
   
    mcmc_data dat;
    coord(ref_path, ss_path, valid_path, \
	  ldmat_path, dat, sz, opt_llk);

    if (dat.beta_mrg.size() == 0) {
	cout << "0 SNPs remained after coordination. Exit." << endl;
	return;
    }
    
    
    ldmat_data ldmat_dat;

    // initialize mcmc
    MCMC_state state = MCMC_state(dat.beta_mrg.size(), M, a0k, b0k, \
	    sz);
    for (size_t i=0; i<dat.beta_mrg.size(); i++) {
	dat.beta_mrg[i] /= c;
    }

    
    MCMC_samples samples = MCMC_samples(dat.beta_mrg.size());
   
    solve_ldmat(dat, ldmat_dat, a, sz, opt_llk);
    state.update_suffstats();

    Function_pool func_pool(n_threads);

    //auto start = steady_clock::now();

    for (int j=1; j<iter+1; j++) {
	state.sample_sigma2();

	for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
	    state.calc_b(i, dat, ldmat_dat);
	}

	for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
	    func_pool.push(std::bind(&MCMC_state::sample_assignment, \
			&state, i, ref(dat), ref(ldmat_dat)));

	    /*func_pool.push(std::bind(&MCMC_state::sample_beta, \
			&state, i, ref(dat), ref(ldmat_dat))); */
	}

	func_pool.waitFinished();

	state.update_suffstats();

	state.sample_V();
	state.update_p();
	state.sample_alpha();

	for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
	    state.sample_beta(i, dat, ldmat_dat);
	}

	//auto start = steady_clock::now();
	state.sample_eta(ldmat_dat);

	if ((j>burn) && (j%thin == 0)) {
	    state.compute_h2(dat);
	    samples.h2 += state.h2*square(state.eta) / n_pst;
	    gsl_blas_daxpy(state.eta/n_pst, state.beta, \
		    samples.beta);
	}

	if (j % 100 == 0) {
	    state.compute_h2(dat);
	    cout << j << " iter. h2: " << state.h2*square(state.eta) << \
		" max beta: " << gsl_vector_max(state.beta)*state.eta \
		 << endl;
		
	}

	//auto stop = steady_clock::now();

	//auto duration = duration_cast<milliseconds>(stop - start);
	//cout << duration.count() << " millisecondss" << endl;
    }

    cout << "h2: " << \
	samples.h2 << " max: " << \
	gsl_vector_max(samples.beta) << endl;

    ofstream out(out_path);
    out << "SNP" << "\t" << \
	"A1" << "\t" << "beta" << endl;

    for (size_t i=0; i<dat.beta_mrg.size(); i++) {
	double tmp = gsl_vector_get(samples.beta, i);
	out << dat.id[i] << "\t" << \
	    dat.A1[i] << "\t" <<  \
	    tmp << endl; 
    }
    out.close();

    for (size_t i=0; i<dat.ref_ld_mat.size(); i++) {
	gsl_matrix_free(ldmat_dat.A[i]);
	gsl_matrix_free(ldmat_dat.B[i]);
	gsl_matrix_free(ldmat_dat.L[i]);
	gsl_vector_free(ldmat_dat.calc_b_tmp[i]);
	gsl_vector_free(ldmat_dat.beta_mrg[i]);
	gsl_matrix_free(dat.ref_ld_mat[i]);
    } 

    return;
}

/*int main() {
    mcmc("test.snpInfo", "/ysm-gpfs/pi/zhao-data//gz222/height/summ_stats/PRS_cs.txt", \
	    "/ysm-gpfs/pi/zhao/gz222/UKB_real/genotype/Ukb_imp_v2_hm3.bim", \
	    "test.dat", "chr1.txt", 252230, 0.1, 1, 1000, \
	    .5, .5, 100, 10, 5, 3);
    return 0;
}*/

