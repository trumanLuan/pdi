
## :::::::::::::::::::::::::::::::::::::::::::::::
## project dirs and data files.
## :::::::::::::::::::::::::::::::::::::::::::::::

rm(list=ls())

proj.dir="E:/projects/PDI"
code.dir <- file.path(proj.dir, "code")
out.dir=file.path(proj.dir, "results")
data.dir = file.path(proj.dir, 'data')
pfm.data.dir <- file.path(data.dir, "pfam/pfam_full_prody")
coevol.data.dir <- file.path(data.dir, "intermediate_dt/coevolution")

## speer.score
load(file.path(data.dir, "intermediate_dt/speer.score.c.rdata") ) # load speer.score.c, which includes speer score of residue sites in MSA


## :::::::::::::::::::::::::::::::::::::::::::::::
## packages and source functions
## :::::::::::::::::::::::::::::::::::::::::::::::

## visualization
require(ggplot2)
require(igraph)

## data and string manipulation
require(stringr)
require(reshape2)
require(dplyr)

## source functions
source( file.path(code.dir, "functions/in_use/readFasta.r") )

## :::::::::::::::::::::::::::::::::::::::::::::::
## overview of dataset for coevolution anlaysis
## 
## :::::::::::::::::::::::::::::::::::::::::::::::
## TF domain sequences used for co-evolution analysis
tf.fam <- c('Homeobox', 'HLH', 'zf-C4', 'Ets', 'Fork_head', 'HMG_box', 'zf-C2H2', 'Zn_clus', 'bZIP_1')
tf.aln.seq.fs <- grep("aln.trim.fasta", dir(pfm.data.dir), value=T)

tf.aln.seq <- foreach(fi = tf.aln.seq.fs, .combine="rbind") %do% {
  tmp.dt <- readFasta(file.path(pfm.data.dir, fi))
  curr.species <- foreach(i = tmp.dt[,1], .combine='c') %do% {
    i.out = stringr::str_split(i, '/')[[1]][1]
    stringr::str_split(i.out, '_')[[1]][2]
    }
  tmp.dt <- data.frame(seq_id=tmp.dt[,1], sequence=tmp.dt[,2], 
                       tf_fam=rep( stringr::str_split(fi, "\\.")[[1]][1], nrow(tmp.dt) ), 
                       species = curr.species, stringsAsFactors = F)
  tmp.dt
}

tf.aln.seq$seq_id <- stringr::str_replace_all(tf.aln.seq$seq_id, "[[:space:]]+", "")

## summary of unique sequences for each tf_fam
tf.aln.seq %>% subset(! species %in% c("0.90", "9ARCH", "9FIRM", "9ORYZ")) %>% group_by(tf_fam) %>% summarize(n_seq = n_distinct(seq_id))

## summary of number of species for each tf fam
tf.aln.seq %>% subset(! species %in% c("0.90", "9ARCH", "9FIRM", "9ORYZ")) %>% group_by(tf_fam) %>% summarize(n_seq = n_distinct(species))




## :::::::::::::::::::::::::::::::::::::::::::::::
## coevolution analysis with the ProDy suite
## four methods were used: MI, MIp, OMES and SCA
## for revised manuscript: DCA, 
## :::::::::::::::::::::::::::::::::::::::::::::::

## ------
## detect coevolving residue pairs using the four different methods.
## run script file
## 1_supp_prody_for_coevolution.py
## linux server, conda env base, ProDy installed in python environment.

## ------
## collect and organize the results from each method

res.fs <- grep("_DCA", dir(coevol.data.dir), value=T)
combine.dt.dca <- foreach(fi = res.fs, .combine="rbind") %do% {
  message(fi)
  tmp.dt <- read.table( file.path(coevol.data.dir, fi), header=F,sep='\t',as.is=T)
  rownames(tmp.dt) = str_c('V',1:nrow(tmp.dt),sep='')
  colnames(tmp.dt) = str_c('V',1:ncol(tmp.dt),sep='')
  
  tmp.dt.melt = melt(as.matrix(tmp.dt )) #melt the matrix derived from a df object.
  sorted.pos = t(apply(tmp.dt.melt, 1, function(x) sort(x[1:2])))
  tmp.dt.melt$Var1 = sorted.pos[,1]
  tmp.dt.melt$Var2 = sorted.pos[,2]
  tmp.dt.melt = unique(tmp.dt.melt)
  tmp.dt.melt = tmp.dt.melt[!(tmp.dt.melt$Var1 == tmp.dt.melt$Var2),]
  
  # the top significant SDS.
  curr.tf.fam <- stringr::str_replace_all(fi, "_DCA.txt", "")
  speer.score = subset(speer.score.c, tf_fam %in% curr.tf.fam)
  speer.score$column = as.integer(speer.score$column)
  speer.score$pval = as.numeric(speer.score$pval)
  sds = str_c('V',subset(speer.score, pval < 0.1)$column+1,sep='')
  
  # co-evolution score for each class: intra-, between-, extra-
  intra.pair = subset(tmp.dt.melt, Var1 %in% sds & Var2 %in% sds)
  between.pair = subset(tmp.dt.melt, (Var1 %in% sds & (!(Var2 %in% sds))) | (Var2 %in% sds & (!(Var1 %in% sds))) )
  extra.pair = subset(tmp.dt.melt, !(Var1 %in% sds | Var2 %in% sds))
  
  # generate data frame used to plot by ggplot.
  intra.pair$class = rep('intra',nrow(intra.pair))
  between.pair$class = rep('between',nrow(between.pair))
  extra.pair$class = rep('extra',nrow(extra.pair))
  
  all.pair = do.call(rbind, list(intra.pair, between.pair, extra.pair))
  all.pair$tf_fam = rep(curr.tf.fam, nrow(all.pair))
  all.pair
}

saveRDS(combine.dt.dca, file=file.path(data.dir, "intermediate_dt/coevol_score_combined_DCA.rds"))


## --------------------
## candidate highly-coevolving-residue pairs (HCR) for each method

## sub-function
re_organize_df <- function(x){ # x is df such as evo.mat.ets
  x$Var1 <- as.numeric(stringr::str_replace(x$Var1, 'V', ''))
  x$Var2 <- as.numeric(stringr::str_replace(x$Var2, 'V', ''))
  tmp <- x[,1:2]
  tmp <- t(apply(tmp, 1, sort))
  x$Var1 <- tmp[,1]; x$Var2 <- tmp[,2]
  x$pairs <- paste(x$Var1, x$Var2, sep='_')
  x
}

## for DCA
combine.dt.dca <- re_organize_df(combine.dt.dca)

combine.dt.dca.top10p <- foreach(tfi = unique(combine.dt.dca$tf_fam), .combine='rbind') %do% {
  message(tfi)
  fi.dt <- combine.dt.dca %>% subset(tf_fam == tfi)
  fi.dt[fi.dt$value > quantile(fi.dt$value, probs = 0.9),]
  
}



## :::::::::::::::::::::::::::::::::::::::::::::::
## for revised manuscript
## compare the defined HCRs to that from the DCA mehtod
## 
## :::::::::::::::::::::::::::::::::::::::::::::::

## //
hcr.df <- read.table(file.path(data.dir, "intermediate_dt/highly_coevolving_reside_pairs.txt"),header=T,sep='\t',as.is=T )

foreach(tfi = unique(hcr.df$tf_fam)) %do% {
  message(tfi)
  tfi.dt <- hcr.df %>% subset(tf_fam == tfi)
  print(table(table(tfi.dt$pairs)) )
  cat("\n")
}

## // how many HCRs obtained using DCA methods can be found in our original manuscript?

foreach(tfi = unique(hcr.df$tf_fam) ) %do% {
  message(tfi)
  
  tfi.hcr.df <- hcr.df %>% subset(tf_fam == tfi)
  tfi.dca <- combine.dt.dca.top10p %>% subset(tf_fam ==tfi)
  
  message("   % of DCA scores in our current results: ", table(tfi.dca$pairs %in% tfi.hcr.df$pairs)['TRUE'] / nrow(tfi.dca) )
  cat("\n")
}

## :::::::::::::::::::::::::::::::::::::::::::::::
## TSDS sites versus highly-coevolving-residue-pairs
## network analysis
## 
## :::::::::::::::::::::::::::::::::::::::::::::::

dt.coevol.score <- hcr.df; # combine.dt.dca.top10p, hcr.df
#~ netowrk demo
i = "HMG_box";
i

tmp <- subset(dt.coevol.score, tf_fam == i, select=c("Var1", "Var2"))
tmp <- unique(tmp)
g <- graph_from_data_frame(tmp, directed=F)

plot(g, edge.arrow.size=.2, edge.curved=0, vertex.color="orange", vertex.frame.color="#555555",
     vertex.label=V(g), vertex.label.color="black", vertex.label.cex=.7, layout=layout_with_fr, main=str_c("Coevolution netowrk: ",i)) # layout_with_lgl, layout_with_fr

## sub graph with degree >2
g.sub <- subgraph(g, V(g)[degree(g)>3]) # Homeobox=4,HLH=3,zf-C4=6,zn_clus=0, Ets=6-7, HMG_box=5-3; bZIP_1=4

# ceb <- cluster_edge_betweenness(g.sub) 
# dendPlot(ceb, mode="hclust")
# plot(ceb, g.sub, main=str_c("coevolved netork: ", i, " hclust")) 

set.seed(123)
cfg <- cluster_fast_greedy(as.undirected(g.sub))
plot(cfg, as.undirected(g.sub), main=str_c("coevolved netork: ", i), layout=layout_with_fr) #, layout=layout_with_fr



