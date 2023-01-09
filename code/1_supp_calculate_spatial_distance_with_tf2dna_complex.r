

#############
## dirs and files
#############

rm(list=ls())
proj.dir="E:/projects/PDI"
code.dir <- file.path(proj.dir, "code")
data.dir = file.path(proj.dir, 'data')
pdb.raw.dir <- file.path(data.dir, "pfam_pdb/pdb")
pdb.aln.dir <- file.path(data.dir, "pfam_pdb/pdb_aln")
pdb.del.dir <- file.path(data.dir, "pfam_pdb/pdb_del")
pbd.mutation.dir <- file.path(data.dir, "pfam_pdb/pdb_mutation")
pdb.better.dir <- file.path(data.dir, "pfam_pdb/pdb_2a_reso")

pdb.meta = read.table(file.path(data.dir,'pfam_pdb/pfam_pdb_all.txt'),header=T,sep='\t',as.is=T)
pdb.meta = subset(pdb.meta, ! TF_family %in% c('zf-C2H2'))
pdb.meta <- pdb.meta %>% subset(PDB_ID %in% stringr::str_replace_all(dir(pdb.better.dir), ".pdb", ""))

aa.alphabet = read.table(file.path(data.dir, 'intermediate_dt/amino_acid_alphabet.txt'), header=T,sep='\t',as.is=T)


#############
## tools and packages
#############

library(Rpdb)
library(stringr)

source(file.path(code.dir, 'functions/in_use/tf2dnaDist.R') )
source(file.path(code.dir, "functions/in_use/readFasta.r") )
source(file.path(code.dir, "functions/in_use/betweenResidueDist.r") )


#########################3
## TF domain between-residue distance calculation.
#########################3

for(curr.tf in unique(pdb.meta$TF_family) ){  
cat(curr.tf, '\n')
curr.pdb.lst = pdb.meta[pdb.meta$TF_family == curr.tf,]

cat('> between-atomes dist for all pdbs...\n')
aa.aa.dist = list()
for(i in 1:nrow(curr.pdb.lst)){
	cat(i,'/',nrow(curr.pdb.lst),'\n')
	x.pdb = read.pdb(file.path(pdb.better.dir, paste(curr.pdb.lst[i,'PDB_ID'],'.pdb',sep='')))
	x.pdb$atoms = x.pdb$atoms[x.pdb$atoms$recname=='ATOM',]
	chain.id = unique(x.pdb$atoms$chainid)
	tf.chain = curr.pdb.lst[i,'PDB_chain_ID']
	if(tf.chain %in% chain.id){ tf2tf.d=betweenResidueDist(x.pdb, tf.chain=tf.chain)	}
	aa.aa.dist[[paste(curr.pdb.lst[i,'PDB_ID'],tf.chain,sep='_')]] = tf2tf.d
} # i end

cat('\n')

##  align protein chains in pdb files to the tf reference.
# extract positions according to the TF reference.
muscle.output = file.path(pdb.better.dir, stringr::str_c('pdb.protein.', curr.tf, ".aln.fa",sep='')) 
pdb.seq.aln = readFasta(muscle.output)
ref.pos = which(strsplit(pdb.seq.aln[pdb.seq.aln[,1]=='TF_REF',2],split='')[[1]] != '-')

# -> calculate index of included postion in each pdb domain.
cat('> Align pdb chains to the reference domain position...\n')
index.vec = NULL
for(i in 1:nrow(pdb.seq.aln)){
	tmp.seq = strsplit(pdb.seq.aln[i,2],split='')[[1]]
	local.index = 0
	for(j in tmp.seq){
		if(j == '-') local.index = c(local.index, 0) # 0 indicates a gap.
		if(j != '-') local.index = c(local.index, max(local.index)+1)
	}
	local.index = local.index[2:length(local.index)]
	keep.index = local.index[ref.pos]
	index.vec = c(index.vec, paste(keep.index, collapse=','))
}

names(index.vec) = pdb.seq.aln[,1]

# -> generate TF-TF distance matrix according to pdb.dist object and index.vec object.
dist.mat = list()
for(i in 1:length(aa.aa.dist)){
	cat(i, '/', length(aa.aa.dist),'\n')
	curr.chain.dist = aa.aa.dist[[i]]
	tmp.label <- labels(aa.aa.dist[i])
	tmp.label <- stringr::str_replace_all(tmp.label, "_", ".")
	keep.index = index.vec[grepl(tmp.label, names(index.vec) ) ]
	
	keep.index = as.integer(strsplit(keep.index, split=',')[[1]])
	local.index = keep.index
	local.index[which(local.index==0)] = length(local.index) + 1
	names(local.index) = 1:length(local.index)

	local.pair = expand.grid(1:length(local.index), 1:length(local.index))
	local.pair = unique(t(apply(local.pair, 1, sort)))
	tmp.d = apply(local.pair, 1,function(x) {
		if(all(c(local.index[x[1]],local.index[x[2]]) < length(local.index))){
			curr.chain.dist[local.index[x[1]],local.index[x[2]]]
			} else {return(NA)}
		})

	dist.mat[[labels(aa.aa.dist[i])]] = data.frame(local.pair, dist=tmp.d)
}

cat('\n')
cat('> mean of between-residues dist...\n')
dist.m = sapply(dist.mat, function(x)x$dist)
dist.m = apply(dist.m, 1, function(x) mean(x,na.rm=T))  # mean of between-residues distance from all PDB entries
dist.m = data.frame(dist.mat[[1]][,1:2],dist=dist.m)
write.table(dist.m, file=file.path(pdb.better.dir, paste(curr.tf, 'between.residues.within.domain.distance',sep='.')), quote=F,sep='\t',col.names=T, row.names=F)
cat('\n')
}


#########################3
## TF-DNA distance calculation.
#########################3

isDNAChain = function(x.pdb, chain.id){
  is.dna = any(x.pdb$atoms[x.pdb$atoms$chainid==chain.id, 'resname'] %in% c('DA','DT','DC','DG')) | 
    any(x.pdb$atoms[x.pdb$atoms$chainid==chain.id, 'resname'] %in% c('A','T','C','G'))
  is.dna
}

# ## pdb TF-DNA interface distance of individual tf domains
for(curr.tf in unique(pdb.meta$TF_family)){ #"Ets" "HMG_box" "Zn_clus" "Homeobox" "HLH" "zf-C4" "bZIP_1"
  cat(curr.tf, '\n')
  curr.pdb.lst = pdb.meta[pdb.meta$TF_family == curr.tf,]
  
  pdb.dist = list() # cal dist between tf residues and DNAs
  for(i in 1:nrow(curr.pdb.lst)){
    cat(i,'in', nrow(curr.pdb.lst),'...\n')
    x.pdb = read.pdb(file.path(pdb.better.dir, paste(curr.pdb.lst[i,'PDB_ID'],'.pdb',sep='')))
    x.pdb$atoms = x.pdb$atoms[x.pdb$atoms$recname=='ATOM',]
    chain.id = unique(x.pdb$atoms$chainid)
    
    is.dna = sapply(chain.id, function(x){ isDNAChain(x.pdb,x) }) # confirm protein chains and DNA chains.
    dna.chain = names(is.dna[is.dna])
    tf.chain = curr.pdb.lst[i,'PDB_chain_ID']
    
    tf2dna.d=sapply(dna.chain, function(x){tf2dnaDist(x.pdb, tf.chain=tf.chain, dna.chain=x)})
    tf2dna.d = apply(tf2dna.d, 1, min)
    
    pdb.dist[[paste(curr.pdb.lst[i,'PDB_ID'],tf.chain,sep='_')]] = tf2dna.d
    
  }
  
  # -> summarize the distance list.
  ifuse = sapply(pdb.dist, length) >0
  pdb.dist = pdb.dist[ifuse]; #length(pdb.dist)
  
  ##  align protein chains in pdb files to the tf reference.
  cat('\n')
  cat('tf-chains-and-reference-seq-alignment...\n')

    # extract positions according to the TF reference.
  muscle.output = file.path(pdb.better.dir, paste('pdb.protein.',curr.tf, ".aln.fa",sep=''))
  pdb.seq.aln = readFasta(muscle.output)
  ref.pos = which(strsplit(pdb.seq.aln[pdb.seq.aln[,1]=='TF_REF',2],split='')[[1]] != '-')
  
  # -> calculate index of included postion in each pdb domain.
  index.vec = NULL
  for(i in 1:nrow(pdb.seq.aln)){
    tmp.seq = strsplit(pdb.seq.aln[i,2],split='')[[1]]
    local.index = 0
    for(j in tmp.seq){
      if(j == '-') local.index = c(local.index, 0) # 0 indicates a gap.
      if(j != '-') local.index = c(local.index, max(local.index)+1)
    }
    local.index = local.index[2:length(local.index)]
    keep.index = local.index[ref.pos]
    index.vec = c(index.vec, paste(keep.index, collapse=','))
  }
  names(index.vec) = pdb.seq.aln[,1]
  
  
  # -> generate TF-DNA distance matrix according to pdb.dist object and index.vec object.
  pdb.dist.vec = unlist(pdb.dist)
  dist.mat = NULL
  for(i in 1:length(pdb.dist)){ # colnames(pdb.lst) = 'use_chain'
    curr.chain.dist = pdb.dist.vec[grepl(labels(pdb.dist[i]), names(pdb.dist.vec))]
    
    tmp.label <- labels(pdb.dist[i])
    tmp.label <- stringr::str_replace(tmp.label, "_", "\\.")
    
    keep.index = index.vec[ grepl(tmp.label, names(index.vec))  ]
    keep.index = as.integer(strsplit(keep.index, split=',')[[1]])
    keep.index[which(keep.index==0)] = length(curr.chain.dist) + 1
    dist.mat = rbind(dist.mat, curr.chain.dist[keep.index])
  }
  
  rownames(dist.mat) = labels(pdb.dist)
  dist.m = colMeans(dist.mat, na.rm=T)
  names(dist.m) = 1:length(dist.m)
  write.table(dist.m, file=file.path(pdb.better.dir, paste(curr.tf, 'tf-dna.distance',sep='.')), quote=F,sep='\t',col.names=F, row.names=T)
  
}

