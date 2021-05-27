
#include <iostream>
#include <string>
#include <vector>
#include "gsl/gsl_matrix.h"

typedef struct {
    std::string A1;
    std::string A2;
    bool include_ref;
    bool include_ss;
    double beta;
    int array;
    double sz;
} CoordInfo;

typedef struct {
    std::vector<std::string> id;
    std::vector<std::string> A1;
    std::vector<std::string> A2;
    std::vector<double> beta_mrg;
    std::vector<std::pair<size_t, size_t>> boundary;
    std::vector<gsl_matrix *> ref_ld_mat;
    std::vector<double> sz;
    std::vector<int> array;
} mcmc_data;

void coord(const std::string &ref_path, const std::string &ss_path, \
	const std::string &valid_path, const std::string &ldmat_path, \
	mcmc_data &dat, unsigned sz, int opt_llk);




