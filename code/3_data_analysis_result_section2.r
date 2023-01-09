
# TF domain sequence MSA, and clustering analysis
##########################
#1. load packages and functions.
rm(list=ls())
setwd('D:/onedrive_cloud/OneDrive/My_Projects/PDI')
data.dir = './data'
out.dir="./results/v6"

source('./code/functions/in_use/readFasta.r')
source('./code/functions/in_use/ToFasta.r')

library(factoextra)
library(NbClust) # for determination of the number of DNA motif clusters
library(mclust)

library(ggplot2)
library(stringr)
library("gridExtra")
library(dplyr)
source("./code/multiplot.r")

##########################
#2. load TF2DNA data.
load(file.path(data.dir, 'TF2DNA/motif_summary_pwm.RData'))
d=d.stamp; tf2dna = sub.sub.tf2dna; rm(d.stamp, sub.sub.tf2dna)
d.matrix = as.matrix(d)

## Description on tf2dna:
## DNA binding motifs - Motif_ID
## TF family name - Family_Name

##########
## 1.TF-binding DNA motifs clustring
## 2. TF SDS prediction using SPEER, Multi-RELIEF, MACHINE-Learning methods such as Random Forest.
## 3. Emsembled scoring of TF residues.
## for all included TF families one by one.
##########

## define motif clusters.
tf.fam = (unique(tf2dna$DBDs))

for(curr.tf.fam in tf.fam){
  # curr.tf.fam = "Zn_clus", Ets, Zn_clus, Fork_head
  cat('#############################\n'); cat(curr.tf.fam, '\n')
  curr.tf.motif = subset(tf2dna, DBDs %in% curr.tf.fam)$Motif_ID
  curr.d.matrix = d.matrix[curr.tf.motif, curr.tf.motif]
  
  p.wss= fviz_nbclust(curr.d.matrix, hcut, k.max=20, method='wss', nboot=50)+labs(subtitle=str_c("Wss stat for ",curr.tf.fam ))
  n.clust = which(p.wss$data$y < p.wss$data$y[1]/10)[1]; # 90% of within-cluster sum of square ratio
  cat(">>> the optimal number of cluster is",n.clust,'\n')
  if(is.na(n.clust)) cat(">>> Large Heterogeneity:", curr.tf.fam, '\n')
  if(is.na(n.clust)) next;
  
  h.clust = hclust(as.dist(curr.d.matrix), method="complete") # get the clusters of hirarchical clustering analysis
  my.clust = cutree(h.clust, k = n.clust)
  cat(">>> the clusters in", curr.tf.fam, ':\n'); cat(table(my.clust))
  cat('\n')
  
  # in.clust.num = names(table(my.clust))[table(my.clust) > 5] # extract the clusters with > 5 members
  # my.clust = my.clust[my.clust %in% in.clust.num]
  my.clust = as.data.frame(my.clust); my.clust$motif.id = rownames(my.clust)
  a = filter(tf2dna, Motif_ID %in% my.clust$motif.id) # command from dplyr
  my.clust$tf.name = a$TF_Name; my.clust$tf.species = a$TF_Species; my.clust$domain.id= a$DBID
  rownames(my.clust) = NULL
  my.clust = arrange(my.clust, my.clust)
  write.table(my.clust, file.path(out.dir, str_c(curr.tf.fam, 'motifs_clusters.txt',sep='-')), quote=F,sep='\t',row.names=F)
  
  ## prepare sorted TF domain MSA fasta file, and clusters size table.
  my.clust.size = as.data.frame(table(my.clust$my.clust))
  my.clust.in.1 = filter(my.clust.size, Freq>4) # include the clusters with >=5 members in SPEER analysis
  # my.clust.in.2 = filter(my.clust.size, Freq<5) # the clusters with <5 members to define a new cluster
  # my.clust = my.clust %>% mutate(my.clust=replace(my.clust, my.clust %in% my.clust.in.2$Var1, 'pool')) %>% as.data.frame()
  my.clust = subset(my.clust, my.clust %in% my.clust.in.1$Var1)
  my.clust$tf2dna.ac = str_c(my.clust$domain.id, my.clust$motif.id,sep='\t')
  my.clust.size = as.data.frame(table(my.clust$my.clust))
  
  tf.msa = readFasta(file.path(data.dir, 'TF2DNA/tf_domain_aln', str_c(curr.tf.fam, '_domain_aln_trim_EXP.fasta', sep='')))
  tf.msa = tf.msa %>% as.data.frame(stringsAsFactors = F)
  
  tf.msa = subset(tf.msa, name.row %in% my.clust$tf2dna.ac)
  tf.msa.sort = tf.msa[order(match(tf.msa$name.row, my.clust$tf2dna.ac)),]
  
  tf.msa.sort.fa = ToFasta(seq.mat = tf.msa.sort[,2], names=tf.msa.sort[,1])
  write.table(tf.msa.sort.fa, file.path(out.dir, str_c(curr.tf.fam, 'domain_MSA_sort.txt',sep='-')), quote=F,sep='\t',row.names=F,col.names=F)
  write.table(my.clust.size$Freq, file.path(out.dir, str_c(curr.tf.fam, 'motifs_matrix_sort_cluster_size.txt',sep='-')),quote=F,sep='\t',row.names=F,col.names=F)
  
}



## 2. TF SDS prediction using SPEER, Multi-RELIEF, for v2.
## SPEER, standalone in winOS
## Multi-RELIEF, website 
##############
# ./SPEER -i cd00423.FSA -o cd00423_bySpeerScore.out -f 1 25 8 
# or ./SPEER -i cd00423.FSA -o cd00423_bySpeerScore.out -pf sizeCluster.file


## 3.compare the results and Emsemble the scores of TF residues.
##############

speer.dir = file.path(out.dir, '1.identification of SDS using tf2dna domain seq/SPEER')
relief.dir = file.path(out.dir, '1.identification of SDS using tf2dna domain seq/Multi-RELIEF')

tf.fam = tf.fam[! tf.fam %in% c('zf-C2H2', 'Fork_head')] # for v2.

speer.score = list()
relief.score = list()

for(curr.tf.fam in tf.fam){
  cat(">>> ", curr.tf.fam, '\n')
  speer.tmp = readLines(file.path(speer.dir, str_c(curr.tf.fam, '_bySpeer_v2.out',sep='')))
  speer.tmp = speer.tmp[which(grepl('#======', speer.tmp))[1]:which(grepl('#======', speer.tmp))[length(which(grepl('#======', speer.tmp)))]]
  speer.content = speer.tmp[5:(length(speer.tmp)-2)]
  speer.content = t(sapply(speer.content, function(x){ str_split(x, '\\ +')[[1]]}))
  rownames(speer.content) = NULL
  colnames(speer.content) = c('column', 'identity', 'score', 'z_score', 'pval')
  speer.score[[curr.tf.fam]] = as.data.frame(speer.content, stringsAsFactors = F)
  
  relief.tmp = read.table(file.path(relief.dir, str_c(curr.tf.fam, '_weights_v2.txt', sep='')), header=T, sep='\t', as.is=T)
  relief.score[[curr.tf.fam]] = relief.tmp$Weight
}

speer.score.c=NULL;
for(i in tf.fam){ a= speer.score[[i]];  a$tf_fam = rep(i, nrow(a)); speer.score.c = rbind(speer.score.c,a)}
relief.score.c=NULL;
for(i in tf.fam){ a= relief.score[[i]]; a = as.data.frame(a); a$tf_fam = rep(i, nrow(a)); relief.score.c = rbind(relief.score.c,a)}
# save(speer.score.c, file=file.path(out.dir,'speer.score.c.Rdata'))
# save(relief.score.c, file=file.path(out.dir, 'relief.score.c.Rdata'))

# ## correlation analysis of scores
# 
# p = list()
# for(i in tf.fam){
#   cat(i, '\n')
#   speer.tmp = subset(speer.score.c, tf_fam%in%i); relief.tmp = subset(relief.score.c, tf_fam%in%i)
#   speer.tmp$column = as.numeric(speer.tmp$column) +1; speer.tmp = arrange(speer.tmp, column)
#   speer.tmp$score = as.numeric(speer.tmp$score);
#   
#   tmp = data.frame(column=speer.tmp$column,speer = speer.tmp$score, relief =relief.tmp[speer.tmp$column,]$a, stringsAsFactors = F )
#   pearson.cor = cor(tmp$speer, tmp$relief, method='pearson')
#   cat('Pearson cor:',pearson.cor, '\n')
#   rank.cor = cor(tmp$speer, tmp$relief, method='spearman')
#   cat('Spearman cor:',rank.cor, '\n')
#   cat('\n>>>>>>>>>>>>>>\n')
#   p[[i]] = ggplot(tmp, aes(x=speer, y=relief))+ geom_point()+ geom_smooth(method=lm)+ggtitle(i)
# }
# 
# grid.arrange(p[[1]],p[[2]],p[[3]], p[[4]], p[[5]], p[[6]],p[[7]], ncol = 3, nrow = 3)


## 4.seq logos of SDS and corresponding motifs
##############

library(motifStack)
tf.msa.dir = 'E:/project/PDI/SPEER_PC/SPEER~1/data-v2-without_pseudo_cluster'

speer=NULL;for(i in names(speer.score)){ tmp=speer.score[[i]];tmp$tf_fam=rep(i, nrow(tmp));speer=rbind(speer, tmp) }
speer$column = as.numeric(speer$column) + 1
speer$pval = as.numeric(speer$pval)

# input tf domain msa data.
for(curr.tf.fam in tf.fam){
  # curr.tf.fam = tf.fam[1]
  cat(curr.tf.fam, 'BEGIN\n')
  curr.msa.data = readFasta(file.path(tf.msa.dir,paste0(curr.tf.fam, '-domain_MSA_sort.txt')))
  clust.size = readLines(file.path(tf.msa.dir, paste0(curr.tf.fam, '-motifs_matrix_sort_cluster_size.txt')));clust.size=as.numeric(clust.size)
  curr.msa = data.frame(tf_name=curr.msa.data[,1], seq=curr.msa.data[,2], sub_family=rep(paste0('sf',1:length(clust.size)),clust.size), stringsAsFactors = F)
  sds.pos = subset(speer, tf_fam%in%curr.tf.fam); sds.pos = subset(sds.pos, pval <0.1)$column ## significant SDS positions
  
  #sds.pos = c(sds.pos, c(30,18,37, 11, 27, 7))
  
  ## motif merge and motifs cluster logos
  plot.motifs = list()
  for(tf.group in unique(curr.msa$sub_family)){ # tf.group is the cluster index.
    cat('SubFamily', tf.group, '\n')
    tmp.motifs = subset(curr.msa, sub_family %in% tf.group)$tf_name
    tmp.motifs = sapply(tmp.motifs, function(x){ str_split(x, '\t')[[1]][2] })  # extract motifs in certain tf clusters
    sub.pfms = pwm.list[tmp.motifs]
    sub.pfms <- lapply(names(sub.pfms), function(.ele, sub.pfms){
      new("pfm",mat=sub.pfms[[.ele]], name=.ele)}, sub.pfms) # create pfm class.
    
    sub.pfms.aln = DNAmotifAlignment(pfms=sub.pfms) # motifs alignment.
    pfmat.aln = lapply(sub.pfms.aln, function(x) x@mat)
    names(pfmat.aln) = sapply(sub.pfms.aln, function(x) x@name)
    pfmat.aln.new = Reduce('+', pfmat.aln)
    pfmat.aln.new = pfmat.aln.new/(length(pfmat.aln))
    motif <- new("pfm", mat=pfmat.aln.new, name=paste(curr.tf.fam,'group',tf.group,length(sub.pfms),'motifs',sep='_'))
    plot.motifs[[paste(curr.tf.fam,'group',tf.group,sep='_')]] = motif
  }
  
  Sys.setenv(R_GSCMD="\"C:/Program Files/gs/gs9.22/bin/gswin64c.exe\"")
  pdf(file.path(out.dir,paste(curr.tf.fam, 'motifs_stack.pdf',sep='_')),width=3,height=1.4*length(plot.motifs))
  plotMotifLogoStack(plot.motifs, ncex=0.8)
  dev.off()
  cat('\n\n')
}

## plot SDS residues for each cluster of TF domains.

for(curr.tf.fam in tf.fam){
  cat(curr.tf.fam, '\n')
  curr.msa.data = readFasta(file.path(tf.msa.dir,paste0(curr.tf.fam, '-domain_MSA_sort.txt')))
  clust.size = readLines(file.path(tf.msa.dir, paste0(curr.tf.fam, '-motifs_matrix_sort_cluster_size.txt')));clust.size=as.numeric(clust.size)
  curr.msa = data.frame(tf_name=curr.msa.data[,1], seq=curr.msa.data[,2], sub_family=rep(paste0('sf',1:length(clust.size)),clust.size), stringsAsFactors = F)
  sds.pos= subset(speer, tf_fam%in%curr.tf.fam); sds.pos = subset(sds.pos, pval <0.1)$column ## significant SDS positions
  
  #sds.pos = c(sds.pos, c(30,18,37, 11, 27, 7))
  
  
  tf.pwm.list =list() # tf domain pwm
  for(tf.group in unique(curr.msa$sub_family)){
    cat(tf.group,'\n')
    tmp.aa.str = readAAStringSet(file.path(tf.msa.dir,paste0(curr.tf.fam, '-domain_MSA_sort.txt')));
    tmp.aa.str = tmp.aa.str[subset(curr.msa, sub_family %in% tf.group)$tf_name,]
    
    tf.pwm = consensusMatrix(tmp.aa.str, as.prob=T) # obtain pwm from tf domain msa.
    tf.pwm = tf.pwm[, sds.pos]; colnames(tf.pwm) = str_c("S",sds.pos,sep='')
    
    tf.pwm.list[[paste(curr.tf.fam,tf.group,length(tmp.aa.str),'motifs',sep='_')]]=tf.pwm
  }
  
  tf.pwm.list <- lapply(names(tf.pwm.list), function(.ele,tf.pwm.list){
    new("pfm",mat=tf.pwm.list[[.ele]], name=.ele, color=colorset(alphabet="AA", colorScheme="chemistry"))}, tf.pwm.list) # c
  pdf(file.path(out.dir,paste(curr.tf.fam, 'and_hCRs_tfs_stack.pdf',sep='_')),width=6,height=2.5*length(tf.pwm.list))
  plotMotifLogoStack(tf.pwm.list, ncex=0.8)
  dev.off()
  
  cat('\n\n')
}

###############/////
## compare between-sds residue pairs frequency to that between-non-sds 
###############/////
## between-residues amino acid pair frequency comparison
## to reveal interdependent relationship between domain residues

library(motifStack)
tf.msa.dir = 'E:/project/PDI/SPEER_PC/SPEER~1/data-v2-without_pseudo_cluster'

freq.df <- NULL
for(curr.tf.fam in tf.fam){
  cat(curr.tf.fam, '\n')
  if(curr.tf.fam %in% c('zf-C2H2', 'Fork_head')) next;
  
  # subfamily data loading.
  curr.msa.data = readFasta(file.path(tf.msa.dir,paste0(curr.tf.fam, '-domain_MSA_sort.txt')))
  clust.size = readLines(file.path(tf.msa.dir, paste0(curr.tf.fam, '-motifs_matrix_sort_cluster_size.txt')));clust.size=as.numeric(clust.size)
  curr.msa = data.frame(tf_name=curr.msa.data[,1], seq=curr.msa.data[,2], sub_family=rep(paste0('sf',1:length(clust.size)),clust.size), stringsAsFactors = F)
  sds.pos = subset(speer, tf_fam%in%curr.tf.fam); 
  sds.pos = subset(sds.pos, pval <0.1)$column ## significant SDS positions
  
  # whole family msa data loading.
  #whole.msa.data <- readFasta(file.path(data.dir, 'TF2DNA/tf_domain_aln', paste0(curr.tf.fam, '_domain_aln_trim_EXP.fasta')))	
  
  # residue pairs to be compared
  msa.seq.mat <- t(sapply(curr.msa$seq, function(x){ str_split(x, '')[[1]]}))
  
  sds.pair <- t(combn(sds.pos, 2))
  n.residue <- nchar(curr.msa$seq[1])
  non.sds.pair <- t(combn(seq(n.residue)[! seq(n.residue) %in% sds.pos],2))
  tmp <- expand.grid(sds.pos, seq(n.residue)[!seq(n.residue) %in% sds.pos])
  non.sds.pair <- rbind(non.sds.pair, as.matrix(tmp)); rm(tmp)
  
  
  for(i in 1:nrow(non.sds.pair)){
    tmp <- paste(msa.seq.mat[,non.sds.pair[i,1]], msa.seq.mat[,non.sds.pair[i,2]])
    tmp.table <- table(tmp)
    tmp <- data.frame(residue.pair = rep(paste0(non.sds.pair[i,1],'vs', non.sds.pair[i,2]), length(tmp.table)),
                      freq = as.matrix(tmp.table)[,1], pair = rownames(as.matrix(tmp.table)),
                      class = rep(c("non.sds"),length(tmp.table)),
                      tf.fam = rep(curr.tf.fam, length(tmp.table)) )
    freq.df <- rbind(freq.df, tmp)
  }
  
  for(i in 1:nrow(sds.pair)){
    tmp <- paste(msa.seq.mat[,sds.pair[i,1]], msa.seq.mat[,sds.pair[i,2]])
    tmp.table <- table(tmp)
    tmp <- data.frame(residue.pair = rep(paste0(sds.pair[i,1],'vs', sds.pair[i,2]), length(tmp.table)),
                      freq = as.matrix(tmp.table)[,1], pair = rownames(as.matrix(tmp.table)),
                      class = rep(c("sds"),length(tmp.table)),
                      tf.fam = rep(curr.tf.fam, length(tmp.table)) )
    freq.df <- rbind(freq.df, tmp)
  }
  #freq.df$tf.fam <- rep(curr.tf.fam, nrow(freq.df))
  
  cat('\n')
}

save(freq.df, file=file.path(out.dir, 'residue_pair_aa_freq.RData'))

##  normalize the freq to ratio for each position-pair
tf.lst <- as.character(unique(freq.df$tf.fam))
freq.df$tf.fam <- as.character(freq.df$tf.fam)
freq.df$residue.pair <- as.character(freq.df$residue.pair)

freq.df.new <- NULL
for(tfi in tf.lst){
  cat(tfi, "...\n")
  tmp <- subset(freq.df, tf.fam %in% tfi)
  
  tmp.pairs <- unique(tmp$residue.pair)
  for(pair.i in tmp.pairs){
    tmp.sub <- subset(tmp, residue.pair %in% pair.i)
    tmp.sub$ratio <- tmp.sub$freq / sum(tmp.sub$freq)
    freq.df.new <- rbind(freq.df.new, tmp.sub)
  }
  
}




# plot freq distribution
library(ggplot2)
ggplot(freq.df.new, aes(x=ratio, colour=class)) + geom_density() + facet_grid(tf.fam ~ .) +
  # scale_x_continuous(limits = c(0, 50)) +
  theme_classic()


ggplot(freq.df.new, aes( x=ratio, colour=class)) + stat_ecdf() + facet_grid(tf.fam ~ .) +
  # scale_x_continuous(limits = c(0, 50)) +
  theme_classic()



###############/////
##// end of "compare between-sds residue pairs frequency to that between-non-sds "
###############/////


###############/////
## CHi-square test of interdependence between SDS 
###############/////
## between-residues amino acid pair frequency comparison
## to reveal interdependent relationship between domain residues

chisquare.res <- NULL
for(curr.tf.fam in tf.fam){
  cat(curr.tf.fam, '\n')
  if(curr.tf.fam %in% c('zf-C2H2', 'Fork_head')) next;
  
  # subfamily data loading.
  curr.msa.data = readFasta(file.path(tf.msa.dir,paste0(curr.tf.fam, '-domain_MSA_sort.txt')))
  clust.size = readLines(file.path(tf.msa.dir, paste0(curr.tf.fam, '-motifs_matrix_sort_cluster_size.txt')));clust.size=as.numeric(clust.size)
  curr.msa = data.frame(tf_name=curr.msa.data[,1], seq=curr.msa.data[,2], sub_family=rep(paste0('sf',1:length(clust.size)),clust.size), stringsAsFactors = F)
  sds.pos = subset(speer, tf_fam%in%curr.tf.fam); 
  sds.pos = subset(sds.pos, pval <0.1)$column ## significant SDS positions
  
  # residue pairs to be compared
  msa.seq.mat <- t(sapply(curr.msa$seq, function(x){ str_split(x, '')[[1]]}))
  
  sds.pair <- t(combn(sds.pos, 2))
  test.df <- NULL
  for(i in 1:nrow(sds.pair)){
    tmp <- msa.seq.mat[,c(sds.pair[i,1], sds.pair[i,2])]
    # tmp <- msa.seq.mat[,1:2] # just for testing
    tmp <- as.data.frame(tmp)	   
    tmp.table <- table(tmp)
    test.res <- chisq.test(tmp.table)
    # fisher.test(as.matrix(tmp.table), simulate.p.value = T)
    test.df <- rbind(test.df, c(residue.pair=paste0(sds.pair[i,1], '_', sds.pair[i,2]), 
                                statistic=test.res$statistic, pval=test.res$p.value,
                                tf.fam=curr.tf.fam))
  }
  chisquare.res <- rbind(chisquare.res, test.df)
  cat('\n')
}

write.table(chisquare.res, file.path(out.dir, 'TF2DNA_SDS_pairs_interdependence_test.txt'),quote=F,sep='\t',row.names=F,col.names=T)


# plot freq distribution
library(ggplot2)
#pdf(file.path(getwd(), out.dir, 'residue_pair_aa_freq.pdf'))
#p = 
ggplot(freq.df, aes(x=freq, colour=class)) + geom_density() + facet_grid(tf.fam ~ .) +
  scale_x_continuous(limits = c(0, 50))
#print(p); 
#dev.off()
###############/////
##// end of "compare between-sds residue pairs frequency to that between-non-sds "
###############/////


## 5. Frequency analysis of combination of tf domain SDSs.
## 
##############

# for residue 53 of homeobox domain MSA.
# for(curr.tf.fam in tf.fam){
curr.tf.fam = 'Homeobox'
cat(curr.tf.fam, '\n')
curr.msa.data = readFasta(file.path(tf.msa.dir,paste0(curr.tf.fam, '-domain_MSA_sort.txt')))
clust.size = readLines(file.path(tf.msa.dir, paste0(curr.tf.fam, '-motifs_matrix_sort_cluster_size.txt')));clust.size=as.numeric(clust.size)
curr.msa = data.frame(tf_name=curr.msa.data[,1], seq=curr.msa.data[,2], sub_family=rep(paste0('sf',1:length(clust.size)),clust.size), stringsAsFactors = F)
sds.pos= subset(speer, tf_fam%in%curr.tf.fam); sds.pos = subset(sds.pos, pval <0.1)$column ## significant SDS positions

freq.matrix = NULL
# 1-order freq distribution in clusters 
seq.matrix = sapply(curr.msa$seq, function(x){ str_split(x, '')[[1]][sds.pos[1]]})
seq.matrix = data.frame(resid = seq.matrix, cluster=curr.msa$sub_family,stringsAsFactors = F)
freq.matrix = table(seq.matrix)

heatmap.2(as.matrix(freq.matrix), Colv=F, Rowv=F, scale='none', 
          col=colorRampPalette(c("white", "red",'black'))(20), key=T, trace='none',sepwidth=c(0.1,0.1), sepcolor="black",
          colsep=1:ncol(as.matrix(freq.matrix)),
          rowsep=1:nrow(as.matrix(freq.matrix)))

# 2-order freq dist in clusters: 46-53
pos.c = unique(expand.grid(sds.pos, sds.pos)) # domain residues combination
pos.c = subset(pos.c, Var1!=Var2)
pos.c = subset(pos.c, Var1 %in% sds.pos[1] | Var2 %in% sds.pos[1])

seq.matrix = sapply(curr.msa$seq, function(x){ str_split(x, '')[[1]][c(pos.c[1,1], pos.c[1,2])]})
seq.matrix = paste(seq.matrix[1,], seq.matrix[2,], sep='')
seq.matrix = data.frame(resid = seq.matrix, cluster=curr.msa$sub_family,stringsAsFactors = F)
freq.matrix.tmp = table(seq.matrix)
# freq.matrix.tmp = rbind()

# 2-order freq dist in clusters: 49-53
seq.matrix = sapply(curr.msa$seq, function(x){ str_split(x, '')[[1]][c(pos.c[2,1], pos.c[2,2])]})
seq.matrix = paste(seq.matrix[1,], seq.matrix[2,], sep='')
seq.matrix = data.frame(resid = seq.matrix, cluster=curr.msa$sub_family,stringsAsFactors = F)
freq.matrix.tmp = table(seq.matrix)

heatmap.2(as.matrix(freq.matrix.tmp), Colv=F, Rowv=F, scale='none', 
          col=colorRampPalette(c("white", "red",'black'))(20), key=T, trace='none',sepwidth=c(0.1,0.1), sepcolor="black",
          colsep=1:ncol(as.matrix(freq.matrix.tmp)),
          rowsep=1:nrow(as.matrix(freq.matrix.tmp)))

# 3-order freq dist in clusters: 53-46-49
seq.matrix = sapply(curr.msa$seq, function(x){ str_split(x, '')[[1]][c(53,46,49)]})
seq.matrix = paste(seq.matrix[1,], seq.matrix[2,], seq.matrix[3,], sep='')
seq.matrix = data.frame(resid = seq.matrix, cluster=curr.msa$sub_family,stringsAsFactors = F)
freq.matrix.tmp = table(seq.matrix)

library(gplots)

heatmap.2(as.matrix(freq.matrix.tmp), Colv=F, Rowv=F, scale='none', 
          col=colorRampPalette(c("white", "red",'black'))(20), key=T, trace='none',sepwidth=c(0.1,0.1), sepcolor="black",
          colsep=1:ncol(as.matrix(freq.matrix.tmp)),
          rowsep=1:nrow(as.matrix(freq.matrix.tmp)))


