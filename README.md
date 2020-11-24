# SDPR

SDPR (Summary statistics based Dirichelt Process Regression) is a method to compute polygenic risk score (PRS) from summary statistics. It is the extension of Dirichlet Process Regression ([DPR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587666/pdf/41467_2017_Article_470.pdf)) to the use of summary statistics. More details can be found in the preprint.

# Installation

We currently test SDPR on linux x86_64. We will check whether the current version works for Mac soon. There are two ways to install SDPR. 

1. download the [precompiled binary](https://github.com/eldronzhou/SDPR/blob/main/bin/SDPR). If you plan to run SDPR on a linux system with a modern intel processor, this version probably works well.

2. compiling from the source. 

```
git clone https://github.com/eldronzhou/SDPR.git
make
```

# Quick start

SDPR can be run from the command line. To see the full list of options, please refer to our manual or type

```bash
SDPR -h
```

SDPR provides two functions: (1) estimating and paritioning of the reference LD matrix (2) perform MCMC to estimate the posterior effect sizes for each SNP. We provide an example usage for the test dataset:

```bash
cd test/

# make the refernce
SDPR -make_ref -ref_prefix genotype/eur_chr22 -chr 22 -ref_dir ref/

# mcmc
SDPR -mcmc -ref_dir ref/ -ss summary_stat/sim_1.txt -N 503 -chr 22 -out result/SDPR_chr22.txt
```

You may refer to the manual for examples of real data applications.

# Help

We provide guidance on how to use SDPR in the manual. If you encounter bugs, request new features, or have any questions related to installing and running SDPR, please report to the [issue](https://github.com/eldronzhou/SDPR/issues) page. 

# License

SDPR is devloped by [Zhao lab](http://zhaocenter.org) at Yale University. The source code is distributed under the [GPL license](https://github.com/eldronzhou/SDPR/blob/main/LICENSE). SDPR uses [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library), which is also redistributed under the GPL license. We provide a copy of [Intel® MKL®](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) library for dynamic linking, which is redistributed under the [Intel Simple Software License](https://github.com/eldronzhou/SDPR/blob/main/MKL/intel-simplified-software-license.pdf).

# Citation

Preprint coming soon.



