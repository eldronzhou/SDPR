
#include <iostream>
#include "LD.h"
#include "mcmc.h"
#include <string.h>
#include "time.h"

using std::cout; using std::endl;
using std::string;

void print_use() {
    cout << "Usage: SDPR -options" << endl << endl
	<< "Example of estimating LD: " << endl 
	<< "	SDPR -make_ref -ref_prefix ./test/1kg -ref_dir ./test/ref -chr 1 -make_ref" << endl << endl
	<< "Example of performing mcmc: " << endl 
	<< "	SDPR -mcmc -ss ./test/ss/ss.txt -ref_dir ./test/ref -chr 1 -out ./test/SDPR_out.txt" << endl << endl
	<< "Full list of options: " << endl << endl
	<< "-make-ref estimate reference LD matrix." << endl << endl
	<< "-mcmc perform MCMC." << endl << endl
	<< " -ss (required for -make_ref) path to the summary statistics file. We recommend to follow our pipeline to clean up summary statistics by running munge_sumstats.py." << endl
	<< " The summary statistics must have the following format (you can change the name of the header line): " << endl << endl  
	<< " SNP	A1	A2	BETA	P" << endl 
	<< " rs737657        A       G       -2.044  0.0409" << endl
	<< " rs7086391       T       C       -2.257  0.024" << endl
	<< " rs1983865       T       C       3.652   0.00026" << endl
	<< " ..." << endl 
	<< " where SNP is the marker name, A1 is the effect allele, A2 is the alternative allele, BETA is the regression coefficient for quantitative traits or log odds ratio for binary traits, and P is the p value." << endl << endl
	<< " -ref_prefix (required for -make_ref) path to the prefix of the bim file for the reference panel, not including the .bim suffix." << endl << endl
	<< " -ref_dir (required) path to the directory that contains the reference LD information output by SDPR, containing .snpInfo and .dat file." << endl << endl
	<< " -valid (optional) path to the bim file for the testing dataset, including the .bim suffix." << endl << endl
	<< " -out (required for -mcmc) path to the output file containing estimated effect sizes." << endl << endl
	<< " -chr (required) chromsome to work on. Currently support 1-22. Recommend to run in parallel." << endl << endl
	<< " -N (required for -mcmc) GWAS sample size." << endl << endl
	<< " -M (optional) Max number of variance components. M must be greater than 4. Default is 1000." << endl << endl
	<< " -opt_llk (optional) Which likelihood to evaluate. 1 for vanilla modified likelihood and 2 for SNPs genotyped on different individuals. Please refer to manuscript or manual for more details. Default is 1." << endl << endl
	<< " -iter (optional) number of iterations for MCMC. Default is 1000." << endl << endl
	<< " -burn (optional) number of burn-in for MCMC. Default is 200." << endl << endl 
	<< " -thin (optional) Thinning for MCMC. Default is 1 (no thin). " << endl << endl
	<< " -n_threads (optional) number of threads to use. Default is 1." << endl << endl
	<< "-r2 (optional) r2 cut-off for parition of independent blocks. Default is 0.1." << endl
	<< " -a (optional) factor to shrink the reference LD matrix. Default is 0.1. Please refer to the manual for more information." << endl << endl
	<< " -c (optional) factor to correct for the deflation. Default is 1. Please refer to the manual for more information." << endl << endl
	<< " -a0k (optional) hyperparameter for inverse gamma distribution. Default is 0.5." << endl << endl
	<< " -b0k (optional) hyperparameter for inverse gamma distribution. Default is 0.5." << endl << endl
	<< " -h print the options." << endl << endl;
}

int main(int argc, char *argv[]) {

    if (argc == 1) {
	print_use();
	return 0;
    }

    string ss_path, ref_prefix, ref_dir, valid_bim, out_path;
    int N = 0, n_threads = 1, chr = 0;
    double a = 0.1, c = 1, a0k = .5, b0k = .5, r2 = .1;
    size_t M = 1000;
    int iter = 1000, burn = 200, thin = 5;
    int make_ref = 0, run_mcmc = 0;
    int opt_llk = 1;

    // pass command line arguments
    int i = 1;
    while (i<argc) {
	char *end;

	if (strcmp(argv[i], "-ss") == 0) {
	    ss_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-ref_prefix") == 0) {
	    ref_prefix = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-ref_dir") == 0) {
	    ref_dir = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-valid") == 0) {
	    valid_bim = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-out") == 0) {
	    out_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-chr") == 0) {   
	    chr = strtol(argv[i+1], &end, 10);
	    if (*end != 0 || chr > 22 || chr < 0) {
		cout << "Incorrect chromosome: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-N") == 0) {
	    N = strtol(argv[i+1], &end, 10);
	    if (*end != '\0' || N <= 0) {
		cout << "Incorrect N: " << argv[i+1] << endl;
		return 0;
	    }
	    if (N <= 1000) {
		cout << "Warning: sample size too small, might" \
		    " not achieve good performance." << endl;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-M") == 0) {
	    M = strtol(argv[i+1], &end, 10);
	    if (*end != 0 || M <= 4) {
		cout << "Incorrect M: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-a") == 0) {
	    a = strtod(argv[i+1], &end);
	    if (*end != '\0' || a < 0) {
		cout << "Incorrect a: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-c") == 0) {
	    c = strtod(argv[i+1], &end);
	    if (*end != '\0' || c <= 0) {
		cout << "Incorrect c: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-a0k") == 0) {
	    a0k = strtod(argv[i+1], &end);
	    if (*end != '\0' || a0k <= 0) {
		cout << "Incorrect a0k: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-b0k") == 0) {
	    b0k = strtod(argv[i+1], &end);
	    if (*end != '\0' || b0k <= 0) {
		cout << "Incorrect b0k: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-iter") == 0) {
	    iter = strtol(argv[i+1], &end, 10);
	    if (*end != '\0' || iter <= 0) {
		cout << "Incorrect iteration: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-burn") == 0) {
	    burn = strtol(argv[i+1], &end, 10);
	    if (*end != '\0' || burn <= 0) {
		cout << "Incorrect number of iterations: " << argv[i+1] << endl;
		return 0;
	    }
	    if (burn >= iter) {
		cout << "Error: burnin is larger than number of iterations." << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-thin") == 0) {
	    thin = strtol(argv[i+1], &end, 10);
	    if (*end != '\0' || thin <= 0) {
		cout << "Incorrect thin: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-n_threads") == 0) {
	    n_threads = strtol(argv[i+1], &end, 10);
	    if (*end != '\0' || n_threads <= 0) {
		cout << "Incorrect number of threads: " << argv[i+1] << endl;
		return 0;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-r2") == 0) {
	    r2 = strtod(argv[i+1], &end);
	    if (r2 <= 0) {
		cout << "Incorrect r2: " << argv[i+1] << endl;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-opt_llk") == 0) {
	    opt_llk = strtol(argv[i+1], &end, 10);
	    if (opt_llk != 1 && opt_llk != 2) {
		cout << "opt_llk must be in 1 or 2." << endl;
	    }
	    i += 2;
	}
	else if (strcmp(argv[i], "-h") == 0) {
	    print_use();
	    return 0;
	}
	else if (strcmp(argv[i], "-make_ref") == 0) {
	    make_ref = 1;
	    i++;
	}
	else if (strcmp(argv[i], "-mcmc") == 0) {
	    run_mcmc = 1;
	    i++;
	}
	else {
	    cout << "Invalid option: " << argv[i] << endl;
	    return 0;
	}
    }

    if (make_ref && run_mcmc) {
	cout << "both -make_ref and -mcmc are specified. Please specify one of them." << endl;
	return 0;
    }

    if (!chr) {
	cout << "Invalid chromosome specified." << endl;
	return 0;
    }

    if (ref_dir.empty()) {
	cout << "Did not specify the directory of reference LD." << endl;
	return 0;
    }

    // compute LD
    if (make_ref) {

	if (ref_prefix.empty()) {
	    cout << "Did not specify the prefix of the bim file for the reference panel." << endl;
	    return 0;
	}

	div_block(ref_prefix, ref_dir, chr, n_threads, r2);
    }

    // mcmc 
    if (run_mcmc) {

	if (ss_path.empty()) {
	    cout << "Did not specify the path to summary statistics." << endl;
	    return 0;
	}

	if (out_path.empty()) {
	    cout << "Did not specify the path of the output file." << endl;
	    return 0;
	}

	if (!N) {
	    cout << "Did not specify GWAS sample size." << endl;
	    return 0;
	}

	string ref_ldmat = ref_dir + "/chr" + \
	       std::to_string(chr) + ".dat";
	string ref_snpinfo = ref_dir + "/chr" + \
	       std::to_string(chr) + ".snpInfo";

	mcmc(ref_snpinfo, ss_path, valid_bim, \
	ref_ldmat, out_path, N, a, c, M, a0k, b0k, iter, \
	burn, thin, n_threads, opt_llk);
    }

    return 0;
}


