#!/usr/bin/bash

# simulate phenotype
gcta64 --bfile ../genotype/eur_chr22 --simu-qt --simu-causal-loci ../effect_size/sim_1.txt --simu-hsq 0.3 --simu-rep 1 --out ./sim_1

# obtain summary statistics
plink2 --bfile ../genotype/eur_chr22 --pheno ./sim_1.phen --glm a0-ref cols=chrom,pos,alt,ref,a1freq,orbeta,se,tz,p,nobs --out ../summary_stat/sim_1

# format conversion
awk '{BEGIN{print "SNP" "A1" "A2" "BETA" "P"} NR>1 {print $3,$6,$4,$9,$12}}' ../summary_stat/sim_1.PHENO1.glm.linear > ../summary_stat/sim_1.txt

# clean-up
rm -f ../summary_stat/sim_1.PHENO1.glm.linear ../summary_stat/*.log ../summary_stat/*.id 
rm -f *.log *.par
