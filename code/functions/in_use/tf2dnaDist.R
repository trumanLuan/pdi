# function to calculate distance of single TF chain to single DNA chain.
tf2dnaDist=function(pdb, tf.chain, dna.chain){
# -> pdb: object of read.pdb().
# -> tf.chain: TF domain chain ID in atoms matrix of pdb.
# -> dna.chain: DNA chain ID in atoms matrix of pdb.
require(Rpdb)
sec1=pdb$atoms$chainid == tf.chain
sec2=pdb$atoms$chainid==dna.chain
d = distances(pdb, sec1, sec2)
d = norm(d, type='xyz')
rownames(d) = pdb$atoms[sec1,'resid']
colnames(d) = pdb$atoms[sec2,'resid']
d.new = matrix(0,nrow=length(unique(rownames(d))),ncol=length(unique(colnames(d))))
rownames(d.new)=unique(rownames(d))
colnames(d.new)=unique(colnames(d))
for(m in 1:length(unique(rownames(d)))){
	m.residue = unique(rownames(d))[m]
	for(n in 1:length(unique(colnames(d)))){
		n.residue = unique(colnames(d))[n]
		tmp.d = d[rownames(d) == m.residue, colnames(d)==n.residue]
		tmp.d = min(tmp.d)
		d.new[m,n]=tmp.d
	}
}
apply(d.new,1,min) # take the minimum distance of TF2DNA as the residue-DNA distance.
}


############################
## determine whether current chain is DNA chain or not.
isDNAChain = function(x.pdb, chain.id){
  is.dna = any(x.pdb$atoms[x.pdb$atoms$chainid==chain.id, 'resname'] %in% c('DA','DT','DC','DG')) | 
    any(x.pdb$atoms[x.pdb$atoms$chainid==chain.id, 'resname'] %in% c('A','T','C','G'))
  is.dna
}


############################
# 

