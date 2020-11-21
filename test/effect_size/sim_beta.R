#!/usr/bin/R

snp = read.table("../genotype/eur_chr22.bim", header=F, stringsAsFactors=F)
snp = snp[,2]
h2s = 0.3; M = length(snp)
group = 100
pi = c(group, M-sum(group)) / M

set.seed(1234)
n_rep = 1
for (h2 in h2s) {
sigma2 = (h2/M) / pi[1]

for (i in 1:n_rep) {
repeat{
id = sample(1:2, size=M, replace=T, prob=pi)
beta = rep(0, M)
beta[id == 1] = rnorm(sum(id == 1), 0, sqrt(sigma2[1]))
if (abs(M*var(beta) - h2) < 1e-3) break
}
res = data.frame(snp, beta)
res = res[res$beta != 0, ]
write.table(res, file=paste0("sim_",i,".txt"), append=F, quote=F, sep="\t", row.names=F,
            col.names=F)
}
}
