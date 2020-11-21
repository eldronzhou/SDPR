
# estimate reference LD
SDPR -make_ref -ref_prefix genotype/eur_chr22 -chr 22 -ref_dir ref/

# mcmc
SDPR -mcmc -ref_dir ref/ -ss summary_stat/sim_1.txt -N 503 -chr 22 -out result/SDPR_chr22.txt
