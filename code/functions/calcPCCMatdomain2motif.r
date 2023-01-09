# correlation of TF and DNA motifs.
# -> PCC was calculated by each AA in TF and each nt in DNA motif

calcPCCMatdomain2motif = function(domain.pwm, motif.pwm, N=10000){
# function to calculate correlation between TF domains and dna motifs.
# -> domain.pwm: tf domain position weight matrix.
# -> pwm.list: a list including all binding motifs of tf indicated in domain.pwm parameter.
# -> N: the number of sampling, used in the process of generation of cross-table.

pcc.mat=c()

for(i in 1:ncol(domain.pwm)){# each AA in tf domain.
	# aa: amino acid 
	aa = domain.pwm[, i] 
	aa.name = names(aa)
	aa.pool = sample(aa.name, N, replace=T, prob=aa)
	pcc.tmp = c()
	
	if(all(dim(motif.pwm) > 0)){# check the validity of motif.pwm, if any(dim(motif.pwm)) == 0
		for(j in 1:ncol(motif.pwm)){ # each nucleotide in tf-binding motif 
		
			# nt: nucleotide
			nt = motif.pwm[, j]
			nt.name = names(nt)
			nt.pool = sample(nt.name, N, replace=T, prob=nt)
			c = table(aa.pool, nt.pool)
			pcc = assocstats(c)$contingency
			pcc.tmp = c(pcc.tmp, pcc)
		}
	}
	
	pcc.mat= cbind(pcc.mat, pcc.tmp)
}
return(pcc.mat)
}

