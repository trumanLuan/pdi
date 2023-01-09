## :::::::::::::::::::::::::::::::::::::::::::::::
## required packages
## :::::::::::::::::::::::::::::::::::::::::::::::

require(dplyr)
require(ggplot2)
library(RColorBrewer)
require(reshape2)

require(MotIV)
require(Rpdb)

## :::::::::::::::::::::::::::::::::::::::::::::::
## project dirs and data files.
## :::::::::::::::::::::::::::::::::::::::::::::::

proj.dir="E:/projects/PDI"
code.dir <- file.path(proj.dir, "code")
out.dir=file.path(proj.dir, "results")
data.dir = file.path(proj.dir, 'data')
tf2dna.data.dir <- file.path(data.dir, "TF2DNA")
tf2dna.tf.seq.dir <- file.path(data.dir, "TF2DNA/tf_domain_aln")
pfam.tf.seq.dir <- file.path(data.dir, "pfam/pfam_full_prody")
pdb.raw.dir <- file.path(data.dir, "pfam_pdb/pdb")
pdb.aln.dir <- file.path(data.dir, "pfam_pdb/pdb_aln")
pdb.del.dir <- file.path(data.dir, "pfam_pdb/pdb_del")
pbd.mutation.dir <- file.path(data.dir, "pfam_pdb/pdb_mutation")
pdb.better.dir <- file.path(data.dir, "pfam_pdb/pdb_2a_reso")

f.tf2dna <- file.path(data.dir, "TF2DNA/tf2dna.experiments.unique.csv")


## ******************************
##DATASET_1: TF2DNA, TF-DNA interaction 

# input raw data
tf2dna=read.csv(f.tf2dna, header=TRUE, as.is=TRUE)

species = c('Homo_sapiens', 'Mus_musculus', 'Drosophila_melanogaster', 'Saccharomyces_cerevisiae', 'Caenorhabditis_elegans' )
tf2dna <- tf2dna %>% subset(TF_Species %in% species)

tf.maj.fam=table(tf2dna[,'DBDs']); 
tf2dna <- tf2dna %>% subset(DBDs %in% names(tf.maj.fam[tf.maj.fam>30]))


## summary of tf2dna dataset: number of tf-dna in each species
table(tf2dna$TF_Species)
tf.species = table(tf2dna$TF_Species)
tf.species = data.frame(species=names(tf.species), count=tf.species, ratio=tf.species/sum(tf.species))
tf.species = tf.species[,c('species','count.Freq','ratio.Freq')]
with(tf.species,pie(count.Freq, labels=paste0(as.character(species), " ", count.Freq), radius=1,
                    col=colorRampPalette(brewer.pal(8,"YlOrBr"))(nrow(tf.species))))


## summary of tf2dna dataset: no. of proteins in each platform
tf.platform = table(tf2dna$Motif_Type)
tf.platform = data.frame(platform=names(tf.platform), count=tf.platform, ratio=tf.platform/sum(tf.platform))
tf.platform = tf.platform[,c('platform','count.Freq','ratio.Freq')]
with(tf.platform,pie(count.Freq, labels=paste0(as.character(platform), " ", count.Freq), radius=1,
                     col=colorRampPalette(brewer.pal(8,"YlOrRd"))(nrow(tf.platform))))


## summary of tf2dna dataset: total number of TFs in each family and the numbers for each species
tf.count <- tf2dna %>% group_by(DBDs, TF_Species) %>% summarise(n=n())
tf.count.matrix <- reshape2::dcast(tf.count, DBDs ~ TF_Species, value.var="n")

tf.count.rowsum <- rowSums(tf.count.matrix[,2:6], na.rm=T)
names(tf.count.rowsum) <- tf.count.matrix$DBDs
tf.count.colsum <- colSums(tf.count.matrix[,2:6], na.rm=T)

tf.count$DBDs <- factor(tf.count$DBDs, levels=names(sort(tf.count.rowsum, decreasing = T)))
tf.count$TF_Species <- factor(tf.count$TF_Species, levels=names(sort(tf.count.colsum, decreasing = T)))
ggplot(data=tf.count, aes(x=DBDs, y=n, fill=TF_Species)) + geom_bar(stat="identity") +
  theme_classic() + 
  scale_fill_manual(values=brewer.pal(5,"Dark2")  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## compare our collected DNA motifs to known database
df.motifs.final <- tf2dna
df.motifs.final <- df.motifs.final[,3:(ncol(df.motifs.final)-1) ]


tf.fam.names <- unique(df.motifs.final$DBDs)
tf.species.names <- unique(df.motifs.final$TF_Species)

motifs.df.overlap.annot <- data.frame(species=NULL, tf_fam=NULL, tf_name=NULL, tf_ID=NULL, annot_motif=NULL)

for(speciesi in tf.species.names){
  cat(speciesi, "...\n")
  df.motifs.final.subset <- subset(df.motifs.final, TF_Species %in% speciesi)
  
  if(! speciesi %in% c("Drosophila_melanogaster", "Caenorhabditis_elegans" )){
    for(tf.j.index in 1:length(df.motifs.final.subset$TF_Name) ){
      tf.j <- df.motifs.final.subset$TF_Name[tf.j.index]
      tmp.motifs <- MotifDb::query(MotifDb, tf.j,  ignore.case=TRUE)
      
      if(length(tmp.motifs)==0) next;
      
      tmp.motifs.annot <- names(tmp.motifs)
      curr.tf.fam <- df.motifs.final.subset[tf.j.index, "DBDs"]
      curr.tf.id <- df.motifs.final.subset[tf.j.index, "DBID"]
      tmp.df <- data.frame(species=rep(speciesi, length(tmp.motifs.annot)), tf_fam=rep(curr.tf.fam, length(tmp.motifs.annot)),
                           tf_name=rep(tf.j, length(tmp.motifs.annot)),
                           tf_ID=rep(curr.tf.id, length(tmp.motifs.annot)), annot_motif=tmp.motifs.annot )
      
      motifs.df.overlap.annot <- rbind(motifs.df.overlap.annot, tmp.df )
      cat(tf.j.index, "in ", nrow(df.motifs.final.subset), "\n")
    }
  }else {
    for(tf.j.index in 1:length(df.motifs.final.subset$DBID) ){
      tf.j <- df.motifs.final.subset$DBID[tf.j.index]
      tmp.motifs <- MotifDb::query(MotifDb, tf.j,  ignore.case=TRUE)
      if(length(tmp.motifs)==0) next;
      
      tmp.motifs.annot <- names(tmp.motifs)
      curr.tf.fam <- df.motifs.final.subset[tf.j.index, "DBDs"]
      curr.tf.name <- df.motifs.final.subset[tf.j.index, "TF_Name"]
      tmp.df <- data.frame(species=rep(speciesi, length(tmp.motifs.annot)), tf_fam=rep(curr.tf.fam, length(tmp.motifs.annot)),
                           tf_name=rep(curr.tf.name, length(tmp.motifs.annot)),
                           tf_ID=rep(tf.j, length(tmp.motifs.annot)), annot_motif=tmp.motifs.annot )
      
      motifs.df.overlap.annot <- rbind(motifs.df.overlap.annot, tmp.df )
      cat(tf.j.index, "in ", nrow(df.motifs.final.subset) , "\n")
    }
  }
}

motifs.df.overlap.annot$annot_DB <- sapply(motifs.df.overlap.annot$annot_motif, function(x){
  stringr::str_split(x, '-')[[1]][2] })

motifs.df.overlap.annot$annot_tf_in_DB <- sapply(motifs.df.overlap.annot$annot_motif, function(x){
  stringr::str_split(x, '-')[[1]][3] })

motifs.df.overlap.annot$annot_tf_in_DB <- sapply(motifs.df.overlap.annot$annot_tf_in_DB, function(x){
  stringr::str_split(x, '_')[[1]][1] })
motifs.df.overlap.annot$annot_tf_in_DB <- sapply(motifs.df.overlap.annot$annot_tf_in_DB, function(x){
  stringr::str_split(x, '\\.')[[1]][1] })

write.table(motifs.df.overlap.annot, file.path(tf2dna.data.dir, "compar_dna_motifs_to_known_motifdb.txt"), quote=F,sep='\t', row.names=F, col.names=T)


## ******************************
## Intermediate DATASET from TF2DNA dataset;

load(file.path(data.dir, "TF2DNA/dna_motifs/motif_summary_pwm.RData"))

d.stamp2=MotIV::motifDistances(pwm.list, DBscores=jaspar.scores, cc="PCC", align="SWU", top=5, go=1, ge=0.5)
# save(list=c("d.stamp","sub.sub.tf2dna","pwm.list"),file=file.path(out.dir,"motif_summary_pwm.RData"))


## ******************************
##DATASET_2: TF domain sequence dataset, corresponding to TF2DNA dataset

source(file.path(code.dir, "functions/readFasta.r") )
source(file.path(code.dir, "functions/ToFasta.r") )
muscle.bin=file.path(code.dir, "tools/muscle3.8.31_i86win32.exe")

dbd=names(table(tf2dna[,'DBDs']))

for(curr.tf in  dbd){
  cat(curr.tf, '\n')
  cat('1.TF2DNA domain sequence multiple alignment...\n')
  tf.fa.file=file.path(tf2dna.tf.seq.dir, paste(curr.tf,"domain_EXP.fa", sep='_'));
  tf.seq = readFasta(tf.fa.file)
  tf.seq[,2] = gsub('\\.','', tf.seq[,2]); 
  tf.seq[,2] = gsub('-','',tf.seq[,2])
  tf.seq[,2] = toupper(tf.seq[,2])
  
  tf.ref.table = read.table(file.path(data.dir, 'tf_domain_position.tsv'),sep='\t', header=T, as.is=T) 
  tf.ref.id = paste(tf.ref.table[tf.ref.table$TF_family == curr.tf, c('TF_ID','Motif_ID')],collapse='\t')
  
  tf.ref.seq = tf.seq[tf.seq[,1] %in% tf.ref.id,]
  
  pfam.fa.aln.file = file.path(tf2dna.tf.seq.dir, paste(curr.tf, '_domain_aln_EXP.fasta',sep=''))
  system(paste(muscle.bin, "-in", tf.fa.file, "-out", pfam.fa.aln.file)) # input file-combined tf ref and pfam seq.
  
  # -> trim aligned pfam sequence according to reference protein.
  # -> keep amino acid positions based on tf.reference.pos.
  cat('2.Keep the ref seq positions\n')
  pfam.fasta.aln = readFasta(pfam.fa.aln.file) # use to extract TF sequence.
  
  # if(tf.ref.id %in% pfam.fasta.aln[,1]){
  pfam.seq.mat = t(sapply(pfam.fasta.aln[,2], function(x)strsplit(x,split='')[[1]])); 
  rownames(pfam.seq.mat)=pfam.fasta.aln[,1]
  keep.pos = which(pfam.seq.mat[tf.ref.id,]!='-')
  pfam.seq.mat.keep = pfam.seq.mat[,keep.pos]
  pfam.seq.keep = apply(pfam.seq.mat.keep, 1, function(x)paste(x,collapse='') )
  if(all(names(pfam.seq.keep)==pfam.fasta.aln[,1])) pfam.fasta.aln[,2]=pfam.seq.keep
  # }
  
  pfam.fasta.trim = ToFasta(seq.mat=pfam.fasta.aln[,2], names=pfam.fasta.aln[,1])
  pfam.fasta.aln.trim.file = file.path(tf2dna.tf.seq.dir, paste(curr.tf, "_domain_aln_trim_EXP.fasta", sep=''))
  write.table(pfam.fasta.trim, pfam.fasta.aln.trim.file, row.names=F, col.names=F, quote=F, sep='\t')
  
  # -> plot familiar logo of TF family.
  cat('3.weblogo...\n')
  weblogo.bin = '/public/software/weblogo/seqlogo'
  system(paste(weblogo.bin, '-f', pfam.fasta.aln.trim.file, '-F PDF -k 0 -p -h 6 -w 30 -n -Y -C 40 -c -o', file.path(tf2dna.tf.seq.dir, paste(curr.tf,'_tf2dna_domain_aln_logo',sep='')), '-t', curr.tf, sep=' '))
  
}


## ******************************
##DATASET_3: a new and independent large-scale TF domain sequences downloaded from pfam database, different from the TF domain sequence from TF2DNA dataset.

source(file.path(code.dir, "functions/readFasta.r") )
source(file.path(code.dir, "functions/ToFasta.r") )
muscle.bin=file.path(code.dir, "tools/muscle3.8.31_i86win32.exe")

dbd=names(table(tf2dna[,'DBDs']))

for(curr.tf in dbd){
  cat(curr.tf, '\n')
  cat('1.Pfam sequence multiple alignment...\n')
  tf.fa.file=file.path(pfam.tf.seq.dir, paste(curr.tf,"pfam_full.fasta", sep='_'));
  tf.seq = readFasta(tf.fa.file)
  tf.seq[,2] = gsub('\\.','', tf.seq[,2]); 
  tf.seq[,2] = gsub('-','',tf.seq[,2])
  tf.seq[,2] = toupper(tf.seq[,2])
  
  tf.ref.file=file.path(tf2dna.tf.seq.dir, paste(curr.tf,"_domain_EXP.fa", sep=''))
  tf.ref.seq = readFasta(tf.ref.file)
  tf.ref.table = read.table( file.path(data.dir, 'tf_domain_position.tsv'),sep='\t', header=T, as.is=T) 
  tf.ref.id = paste(tf.ref.table[tf.ref.table$TF_family == curr.tf, c('TF_ID','Motif_ID')],collapse='\t')
  tf.ref.seq = tf.ref.seq[tf.ref.seq[,1] %in% tf.ref.id,]
  
  #pfam.seq = rbind(tf.seq, pfam.seq) ###### here
  pfam.seq = rbind(tf.ref.seq, tf.seq)
  pfam.seq.fa = ToFasta(seq.mat=pfam.seq[,2], names=pfam.seq[,1])
  pfam.seq.fa.new = file.path(pfam.tf.seq.dir, paste(curr.tf, ".ualn.new.fasta", sep=''))
  write.table(pfam.seq.fa, pfam.seq.fa.new, row.names=F, col.names=F, quote=F, sep='\t')
  
  pfam.fa.aln.file = file.path(pfam.tf.seq.dir, paste(curr.tf, '.aln.fasta',sep=''))
  system(paste(muscle.bin, "-in", pfam.seq.fa.new, "-out", pfam.fa.aln.file)) # input file-combined tf ref and pfam seq.
  
  # -> trim aligned pfam sequence according to reference protein.
  # -> keep amino acid positions based on tf.reference.pos.
  cat('2.align to pfam...\n')
  pfam.fasta.aln = readFasta(pfam.fa.aln.file) # use to extract TF sequence.
  
  if(tf.ref.id %in% pfam.fasta.aln[,1]){
    pfam.seq.mat = t(sapply(pfam.fasta.aln[,2], function(x)strsplit(x,split='')[[1]])); 
    rownames(pfam.seq.mat)=pfam.fasta.aln[,1]
    keep.pos = which(pfam.seq.mat[tf.ref.id,]!='-')
    pfam.seq.mat.keep = pfam.seq.mat[,keep.pos]
    pfam.seq.keep = apply(pfam.seq.mat.keep, 1, function(x)paste(x,collapse='') )
    if(all(names(pfam.seq.keep)==pfam.fasta.aln[,1])) pfam.fasta.aln[,2]=pfam.seq.keep
  }
  
  pfam.fasta.trim = ToFasta(seq.mat=pfam.fasta.aln[,2], names=pfam.fasta.aln[,1])
  pfam.fasta.aln.trim.file = file.path(pfam.tf.seq.dir, paste(curr.tf, ".aln.trim.fasta", sep=''))
  write.table(pfam.fasta.trim, pfam.fasta.aln.trim.file, row.names=F, col.names=F, quote=F, sep='\t')
  
  # -> plot familiar logo of TF family.
  cat('3.weblogo...\n')
  weblogo.bin = '/public/software/weblogo/seqlogo'
  system(paste(weblogo.bin, '-f', pfam.fasta.aln.trim.file, '-F PDF -k 0 -p -h 6 -w 30 -n -Y -C 40 -c -o', file.path(pfam.tf.seq.dir,  paste(curr.tf,'_seq_logo_whole_family',sep='')), '-t', curr.tf, sep=' '))
}


## ******************************
##DATASET_4: TF-DNA strucutre dataset, downloaded from the PDB database.

## downlaod structure data: .pdb files 
pfam.pdb = read.table( file.path(pdb.raw.dir, '../pfam_pdb_all.txt'), header=T,sep='\t', as.is=T)
pdb.lst = unique(pfam.pdb$PDB_ID)

# for(pdb.i in pdb.lst){
#   message(pdb.i, " ", which(pdb.lst == pdb.i))
#   download.file(url=paste('http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=',pdb.i,sep=''),
#                 destfile=file.path('E:/', paste(pdb.i,'.pdb',sep='')))
# }

# remove pdb files corresponding to those excluded TF families.
rm.pdb.lst = unique(subset(pfam.pdb, TF_family %in% 'zf-C2H2')$PDB_ID)

for(i in dir(pdb.raw.dir)){
  cat(i, '\n')
  if(i %in% str_c(rm.pdb.lst, '.pdb', sep='')){
    file.copy(from=file.path(pdb.raw.dir, i), to=file.path(pdb.del.dir, i))
    file.remove(file.path(pdb.raw.dir,i)) }
}


## remove pdbs which contain only DNA chains or protein chains(mainly).
isDNAChain = function(x.pdb, chain.id){
  is.dna = any(x.pdb$atoms[x.pdb$atoms$chainid==chain.id, 'resname'] %in% c('DA','DT','DC','DG')) | 
    any(x.pdb$atoms[x.pdb$atoms$chainid==chain.id, 'resname'] %in% c('A','T','C','G'))
  is.dna
}

for(i in dir(pdb.raw.dir)){
  cat(i,'\n')
  tmp = readLines(file.path(pdb.raw.dir, i))
  end.index = which(grepl('ENDMDL     ', tmp))
  if(length(end.index)<1){
    
  }else{
    tmp = tmp[1:end.index[length(end.index)]]
    write.table(tmp, file.path(pdb.raw.dir, i),row.names=F, col.names=F, quote=F,sep='\t')
  }
  
  x.pdb = Rpdb::read.pdb(file.path(pdb.raw.dir, i))
  x.pdb$atoms = x.pdb$atoms[x.pdb$atoms$recname=='ATOM',]
  chain.id = unique(x.pdb$atoms$chainid)
  if(!any(sapply(chain.id, function(x) {isDNAChain(x.pdb, x)}))){
    # copy to pdb.del.dir
    file.copy(from=file.path(pdb.raw.dir, i), to=file.path(pdb.del.dir, i))
    file.remove(file.path(pdb.raw.dir,i))
  } 
  cat(which(dir(pdb.raw.dir) == i),'\n')
}

## ...................................................
##  include pdbs, of which the resolution < 4 A.
## remove pdb structures with the resulution > 4A and without resolution information.
for(i in dir(pdb.raw.dir)){
  cat(i,'\n')
  x.pdb = Rpdb::read.pdb(file.path(pdb.raw.dir, i))
  x.pdb.remark = x.pdb$remark
  x.pdb.resolution = grep("REMARK   2 RESOLUTION", x.pdb.remark, value=T)
  x.pdb.resolution2 = stringr::str_extract(x.pdb.resolution, "\\d+\\.\\d+ ANGSTROMS")
  x.pdb.resolution2 = as.numeric(stringr::str_extract(x.pdb.resolution2, "\\d+\\.\\d+"))
  
  if(is.na(x.pdb.resolution2)){
    file.copy(from=file.path(pdb.raw.dir, i), to=file.path(pdb.del.dir, i))
    file.remove(file.path(pdb.raw.dir,i));
    next
  }
  
  if(x.pdb.resolution2 > 4){
    file.copy(from=file.path(pdb.raw.dir, i), to=file.path(pdb.del.dir, i))
    file.remove(file.path(pdb.raw.dir,i))
  }
  cat(which(dir(pdb.raw.dir) == i),'\n')
}

## prepare .fasta files according to the chains in pdb strucutres
source(file.path(code.dir, "functions/generateSeqFromPDB.r") )
source(file.path(code.dir, "functions/readFasta.r") )
source(file.path(code.dir, "functions/ToFasta.r") )

aa.alphabet = read.table(file.path(data.dir, 'amino_acid_table.txt'), header=T,sep='\t',as.is=T)

pdb.files = dir(pdb.raw.dir)
pdb.lst = pdb.lst[ stringr::str_c(pdb.lst,'.pdb',sep='')%in%pdb.files ]

pdb.table <- read.table(file.path(data.dir, "pfam_pdb/pfam_pdb_all.txt"), header=T, sep='\t', as.is=T)
pdb.table <- subset(pdb.table, PDB_ID %in% pdb.lst)

pdb.prot.fa <- NULL
for(pdbi in unique( pdb.lst )){
  cat(pdbi, which(unique(pdb.lst) == pdbi ), "\n")
  tmp <- subset(pdb.table, PDB_ID %in% pdbi)
  
  for(j in 1:nrow(tmp)){
    tmp.sub <- tmp[j,]
    tmp.pdb <- 	Rpdb::read.pdb(file.path(pdb.raw.dir, paste(tmp.sub$PDB_ID,'.pdb',sep='')))
    tmp.chain <- tmp.sub$PDB_chain_ID
    tmp.fa <- generateSeqFromPDB(tmp.pdb, chain.id=tmp.chain, aa.alphabet=aa.alphabet)
    pdb.prot.fa <- rbind(pdb.prot.fa, data.frame(pdb_chain=stringr::str_c(tmp.sub$TF_family, ".", tmp.sub$UniProt_entry, ".", tmp.sub$PDB_ID, ".", tmp.sub$PDB_chain_ID, ".", 
                                                                          stringr::str_replace_all(tmp.sub$PDB_residues, ' ', '') ), seq=tmp.fa, stringsAsFactors = F))
  }
}

## split .fasta files for each TF family
pdb.prot.fa$TF_fam <- sapply(pdb.prot.fa$pdb_chain, function(x) stringr::str_split(x, "\\.")[[1]][1])
  
for(tfi in unique(pdb.prot.fa$TF_fam)){
  tmp <- subset(pdb.prot.fa, TF_fam %in% tfi)
  fa.from.tmp <- ToFasta(tmp$seq, names = tmp$pdb_chain)
  write.table(fa.from.tmp, file.path(pdb.aln.dir, stringr::str_c("pdb.protein.", tfi, ".fa")), quote=F,sep='\t', row.names = F, col.names = F)
}

## mutiple alignment using MUSCLE
muscle.bin=file.path(code.dir, "tools/muscle3.8.31_i86win32.exe")

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.Ets.aln.fa") # for Ets
fa.input = file.path(pdb.aln.dir, "pdb.protein.Ets.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.HMG_box.aln.fa") # for HMG_box
fa.input = file.path(pdb.aln.dir, "pdb.protein.HMG_box.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.Zn_clus.aln.fa") # for Zn_clus
fa.input = file.path(pdb.aln.dir, "pdb.protein.Zn_clus.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.Homeobox.aln.fa") # for Homeobox
fa.input = file.path(pdb.aln.dir, "pdb.protein.Homeobox.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.HLH.aln.fa") # for HLH
fa.input = file.path(pdb.aln.dir, "pdb.protein.HLH.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.zf-C4.aln.fa") # for zf_c4
fa.input = file.path(pdb.aln.dir, "pdb.protein.zf-C4.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.Homeobox.aln.fa") # for Homeobox
fa.input = file.path(pdb.aln.dir, "pdb.protein.Homeobox.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.aln.dir, "pdb.protein.bZIP_1.aln.fa") # for bZIP_1
fa.input = file.path(pdb.aln.dir, "pdb.protein.bZIP_1.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

## ...................................................
##  include pdbs, of which the resolution <= 2.5A .
pdb.meta <- pfam.pdb
pdb.meta <- pdb.meta %>% subset(PDB_ID %in% stringr::str_replace_all( dir(pdb.raw.dir), "\\.pdb", '' ))
pdb.meta$resolution <- rep(NA, nrow(pdb.meta))

for(i in dir(pdb.raw.dir)){
  cat(i,'\n')
  x.pdb = Rpdb::read.pdb(file.path(pdb.raw.dir, i))
  x.pdb.remark = x.pdb$remark
  x.pdb.resolution = grep("REMARK   2 RESOLUTION", x.pdb.remark, value=T)
  x.pdb.resolution2 = stringr::str_extract(x.pdb.resolution, "\\d+\\.\\d+ ANGSTROMS")
  x.pdb.resolution2 = as.numeric(stringr::str_extract(x.pdb.resolution2, "\\d+\\.\\d+"))
  
  pdb.meta$resolution[which(pdb.meta$PDB_ID == stringr::str_replace_all(i, "\\.pdb", ''), arr.ind=T)] <- x.pdb.resolution2
}

pdb.meta.better <- pdb.meta %>% subset(resolution <= 2.5)

## prepare .fasta files according to the chains in pdb strucutres
source(file.path(code.dir, "functions/in_use/generateSeqFromPDB.r") )
source(file.path(code.dir, "functions/in_use/readFasta.r") )
source(file.path(code.dir, "functions/in_use/ToFasta.r") )

aa.alphabet = read.table(file.path(data.dir, 'intermediate_dt/amino_acid_alphabet.txt'), header=T,sep='\t',as.is=T)

# pdb.to.copy <- unique(pdb.meta.better$PDB_ID)
# sapply(pdb.to.copy, function(f) {
#   f.name <- stringr::str_c(f, ".pdb")
#   file.copy(from =file.path(pdb.raw.dir, f.name), to = file.path(pdb.better.dir, f.name))
# })

pdb.files = dir(pdb.better.dir)
pdb.lst = unique(pdb.meta.better$PDB_ID)

# pdb.table <- read.table(file.path(data.dir, "pfam_pdb/pfam_pdb_all.txt"), header=T, sep='\t', as.is=T)
# pdb.table <- subset(pdb.table, PDB_ID %in% pdb.lst)

pdb.prot.fa <- NULL
for(pdbi in unique( pdb.lst )){
  cat(pdbi, which(unique(pdb.lst) == pdbi ), "\n")
  tmp <- subset(pdb.meta.better, PDB_ID %in% pdbi)
  
  for(j in 1:nrow(tmp)){
    tmp.sub <- tmp[j,]
    tmp.pdb <- 	Rpdb::read.pdb(file.path(pdb.better.dir, paste(tmp.sub$PDB_ID,'.pdb',sep='')))
    tmp.chain <- tmp.sub$PDB_chain_ID
    tmp.fa <- generateSeqFromPDB(tmp.pdb, chain.id=tmp.chain, aa.alphabet=aa.alphabet)
    pdb.prot.fa <- rbind(pdb.prot.fa, data.frame(pdb_chain=stringr::str_c(tmp.sub$TF_family, ".", tmp.sub$UniProt_entry, ".", tmp.sub$PDB_ID, ".", tmp.sub$PDB_chain_ID, ".", 
                                                                          stringr::str_replace_all(tmp.sub$PDB_residues, ' ', '') ), seq=tmp.fa, stringsAsFactors = F))
  }
}

## split fasta according to TF family
## add TF domain reference sequence into fasta files.
pdb.prot.fa$tf_fam <- sapply(pdb.prot.fa$pdb_chain, function(x) stringr::str_split(x, "\\.")[[1]][1])

for(tfi in unique(pdb.prot.fa$tf_fam)){
  message(tfi)
  tmp <- subset(pdb.prot.fa, tf_fam %in% tfi)
  tmp.old.fa <- readFasta(file.path(pdb.aln.dir, stringr::str_c("pdb.protein.", tfi, ".fa")))
  tmp.old.fa <- tmp.old.fa[tmp.old.fa[,1] == 'TF_REF',]
  tmp <- rbind(data.frame(pdb_chain = "TF_REF", seq=tmp.old.fa[2], tf_fam = tfi, stringsAsFactors = F), tmp)
  fa.from.tmp <- ToFasta(tmp$seq, names = tmp$pdb_chain)
  write.table(fa.from.tmp, file.path(pdb.better.dir, stringr::str_c("pdb.protein.", tfi, ".fa")), quote=F,sep='\t', row.names = F, col.names = F)
}


## mutiple alignment using MUSCLE
muscle.bin=file.path(code.dir, "tools/muscle3.8.31_i86win32.exe")

fa.aln =  file.path(pdb.better.dir, "pdb.protein.Ets.aln.fa") # for Ets
fa.input = file.path(pdb.better.dir, "pdb.protein.Ets.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.HMG_box.aln.fa") # for HMG_box
fa.input = file.path(pdb.better.dir, "pdb.protein.HMG_box.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.Zn_clus.aln.fa") # for Zn_clus
fa.input = file.path(pdb.better.dir, "pdb.protein.Zn_clus.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.Homeobox.aln.fa") # for Homeobox
fa.input = file.path(pdb.better.dir, "pdb.protein.Homeobox.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.HLH.aln.fa") # for HLH
fa.input = file.path(pdb.better.dir, "pdb.protein.HLH.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.zf-C4.aln.fa") # for zf_c4
fa.input = file.path(pdb.better.dir, "pdb.protein.zf-C4.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.Homeobox.aln.fa") # for Homeobox
fa.input = file.path(pdb.better.dir, "pdb.protein.Homeobox.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 

fa.aln =  file.path(pdb.better.dir, "pdb.protein.bZIP_1.aln.fa") # for bZIP_1
fa.input = file.path(pdb.better.dir, "pdb.protein.bZIP_1.fa")
system(paste(muscle.bin, "-in", fa.input, "-out", fa.aln)) 


## calculate TF domain between-residue 3D distance with the Rpdb package
## // see code in 1_supp_tf2tf_residue_distance.r


## calculate TF-DNA spatial distance with the Rpdb package.
## // see code in 1_supp_tf2dna_distance.r


