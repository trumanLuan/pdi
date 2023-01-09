# code to extract protein chain seq from pdb file.
generateSeqFromPDB = function(x.pdb, chain.id, aa.alphabet=aa.alphabet){
	chain.ids = unique(x.pdb$atoms$chainid)
	if(! chain.id %in% chain.ids) cat('Your defined chainid is not contained in the structure!\n')
	if(chain.id %in% chain.ids){
		sub.atoms = x.pdb$atoms[x.pdb$atoms$chainid %in% chain.id & x.pdb$atoms$recname %in% "ATOM" ,]
		
		residues = unique(sub.atoms$resid)
		chain.seq = sapply(residues, function(x) {unique(sub.atoms[sub.atoms$resid==x,'resname'])} )
		chain.seq = sapply(chain.seq, function(x){aa.alphabet[aa.alphabet[,'long_name']%in%x,'short_name']})
	}
	paste(chain.seq,collapse='')
}
