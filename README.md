
# SDPR

SDPR (Summary statistics based Dirichelt Process Regression) is a method to compute polygenic risk score (PRS) from summary statistics. It is the extension of Dirichlet Process Regression ([DPR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587666/pdf/41467_2017_Article_470.pdf)) to the use of summary statistics. More details can be found in the preprint.


## Installation

Currently we have only tested SDPR on linux x86 as SDPR explictly uses SSE instruction. To install SDPR, you need to first download the repo:

```
git clone https://github.com/eldronzhou/SDPR.git
```

There are two ways to install SDPR, assuming you are working on linux x86. 

1. If you plan to run SDPR on a linux system with a modern intel processor, you may use the precompiled binary `SDPR`. Please make sure that dynamic libraries `gsl/lib/libgsl.so` and `MKL/lib/libmkl_rt.so` are not changed, otherwise SDPR is unable to load the libraries. If you are not able to run the MKL library, you can use openBLAS instead as described in the section below. If you encounter errors reporting cannot load these two dynamic libraries, you may want to try `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:gsl/lib/:MKL/lib/`.

2. If you want to compile SDPR from the source for best performance, you need a C/C++ compiler like g++ (tested under version 4.8.5), GSL (version 2.60) and MKL library (version 2017). For convenience, we redistribute the compiled GSL and MKL building on our Intel Xeon Processors in the repo. To install, type `make`. If this version does not work, please report the error to the issue page. If the issue is related to GSL, you may want to download the source code of GSL and compile it yourself. For details about downloading and installing GSL, please refer to [this page](https://www.gnu.org/software/gsl/) and [this page](https://www.gnu.org/software/gsl/doc/html/usage.html#compiling-and-linking). A tutorial of installing gsl is given on this [page](https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/). If you have other versions of MKL library, please refer to [this manual](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html) for linking advice. On AMD CPU, we recommend using OpenBLAS, please see [here](https://github.com/eldronzhou/SDPR/issues/3) for instructions.

## Quick start

SDPR can be run from the command line. To see the full list of options, please type

```bash
./SDPR -h
```

SDPR provides two functions: (1) estimating and paritioning the reference LD matrix (2) perform MCMC to estimate the posterior effect sizes for each SNP. We provide an example usage for the test dataset:

```bash
cd test/

# make the refernce
../SDPR -make_ref -ref_prefix genotype/eur_chr22 -chr 22 -ref_dir ref/

# mcmc
../SDPR -mcmc -ref_dir ref/ -ss summary_stat/sim_1.txt -N 503 -chr 22 -out result/SDPR_chr22.txt
```

## Input 

### Reference LD

If you are working on summary statistics from EUR ancestry, you can download the [MHC included](https://drive.google.com/file/d/1BhZdUNG9GfCaeH3m0Ou7yE2Ljvp39-j7/view?usp=sharing) or [MHC removed](https://drive.google.com/file/d/1o-cqLDHfqouribbzlEjfM0ikedQ9H4jc/view?usp=sharing) reference LD directory consisting 1 million HapMap3 SNPs estimated from 503 1000 Genome EUR samples. You can also create the reference LD for your preferred reference panel by running the command below. 

```
# 1 thread 
for i in {1..22}; do
./SDPR -make_ref -ref_prefix prefix_ref_genotype -chr ${i} -ref_dir SDPR_ref 
done

# parallel over chr and using 3 threads for each chromosome is recommended
```

If you have memory issues regarding the large size of blocks, you can increase the value of `-r2` to 0.3 or 0.4. This will make the size of independent blocks smaller. 

### Summary Statistics 

The summary statistics should at least contain following columns with the same name (order of the column is not important).

```
SNP	A1	A2	BETA	P
rs737657        A       G       -2.044  0.0409
rs7086391       T       C       -2.257  0.024
rs1983865       T       C       3.652   0.00026
...
```

where SNP is the marker name, A1 is the effect allele, A2 is the alternative allele, BETA is the regression coefficient for quantitative traits or log odds ratio for binary traits, and P is the p value. 

The users may also wish to provide Z scores instead of P value and BETA for input. This is recommended if P values for some variants are too small.

```
SNP	A1	A2	Z
rs737657        A       G       -2.044  
rs7086391       T       C       -2.257  
rs1983865       T       C       3.652   
...
```

We recommend to include the column N if the sample size for each SNP is available, for example

```
SNP     A1      A2      BETA    P       N
rs737657        A       G       -2.044  0.0409      252156
rs7086391       T       C       -2.257  0.024      248425
rs1983865       T       C       3.652   0.00026    253135
...
```

Please refer to the [manual](http://htmlpreview.github.io/?https://github.com/eldronzhou/SDPR/blob/main/doc/Manual.html) for preprocessing of public GWAS summary statistics.

## Running SDPR

When running in parallel using (22*3 = 66 threads), SDPR is able to finish MCMC in around 15 minutes for ~0.8 million Hapmap3 SNPs on an Intel Xeon processor.

```
# for chr i using 3 threads
# recommend to run each chr in parallel
./SDPR -mcmc -ref_dir ref/ -ss ss.txt -N 10000 -chr i -out res_i.txt -n_threads 3
```

The output has format:


```
SNP     A1      beta
rs12255619      C       -0.000124535
rs7909677       G       -0.000106013
rs10904494      C       -0.000178207
...
```

where SNP is the marker name, A1 is the effect allele, beta is the estimated posterior effect sizes.

Once having the ouput, one can use [PLINK](https://www.cog-genomics.org/plink/1.9/score) to derive the PRS.

```
plink --bfile test_geno --score res.txt 1 2 3 header --out test
```

## Help

We provide detailed guidance on how to use SDPR in the [manual](http://htmlpreview.github.io/?https://github.com/eldronzhou/SDPR/blob/main/doc/Manual.html). If you encounter bugs, request new features, or have any questions related to installing and running SDPR, please report to the [issue](https://github.com/eldronzhou/SDPR/issues) page. 

## License

SDPR is devloped by [Zhao lab](http://zhaocenter.org) at Yale University. The source code is distributed under the [GPL license](https://github.com/eldronzhou/SDPR/blob/main/LICENSE). SDPR uses [GSL](https://www.gnu.org/software/gsl/) (GNU Scientific Library), which is also redistributed under the GPL license. We provide a copy of [Intel® MKL®](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) library for dynamic linking, which is redistributed under the [Intel Simple Software License](https://github.com/eldronzhou/SDPR/blob/main/MKL/intel-simplified-software-license.pdf).

## Citation

Zhou G, Zhao H. A fast and robust Bayesian nonparametric method for prediction of complex traits using summary statistics. PLoS Genet. 2021 Jul 26;17(7):e1009697. doi: 10.1371/journal.pgen.1009697. 



