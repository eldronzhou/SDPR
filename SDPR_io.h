
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_matrix.h"

typedef struct {
    std::vector<std::string> id;
    std::vector<std::string> A1;
    std::vector<std::string> A2;
    std::vector<size_t> Pos;
    std::map<size_t, size_t> chr_idx;
} SnpInfo;

size_t get_nsamples(const std::string &fam_path);

void read_bim(const std::string &bim_path, SnpInfo *snpinfo);

void read_bed(gsl_matrix *snp, const std::string &bed_path, const size_t n_samples, size_t left, size_t right);
