## :::::::::::::::::::::::::::::::::::::::::::::::
## project dirs and data files.
## :::::::::::::::::::::::::::::::::::::::::::::::

rm(list=ls())

proj.dir="D:/onedrive_cloud/OneDrive/1My_Documents/projects/my_PDI"
code.dir <- file.path(proj.dir, "code")
out.dir=file.path(proj.dir, "results")
data.dir = file.path(proj.dir, 'data')
pfm.data.dir <- file.path(data.dir, "pfam/pfam_full_prody")
coevol.data.dir <- file.path(data.dir, "intermediate_dt/coevolution")

## speer.score
load(file.path(data.dir, "intermediate_dt/speer.score.c.rdata") ) # load speer.score.c, which includes speer score of residue sites in MSA
library(foreach)
library(ggplot2)
library(ggpubr)

## :::::::::::::::::::::::::::::::::::::::::::::::
## load data.
## :::::::::::::::::::::::::::::::::::::::::::::::

##
prody.root <- coevol.data.dir
mistic.root <- file.path(data.dir, "intermediate_dt/mistic_v1_results")

tf.fams <- c("bZIP_1", 'Ets', 'Fork_head', 'HLH', 'HMG_box', 'Homeobox', 'zf-C2H2', 'zf-C4', 'Zn_clus')

mi.list <- list()

for(tfi in tf.fams) {
  print(tfi)
  
  print('__Input MISTIC MI results:')
  tfi.mi <- read.table(file.path(mistic.root, stringr::str_c(tfi,"_mistic/MI_data")), sep='\t', header=F, as.is=T)
  colnames(tfi.mi) <- c("resnum1", 'resname1', 'resnum2', 'resname2', "mistic_value")
  
  print('__Input prody MI results:')
  tfi.mi.another <- read.table(file.path(prody.root, stringr::str_c(tfi,"_MI.txt")), sep='\t', header=F, as.is=T)
  tfi.mi$prody_value <- foreach(i=1:nrow(tfi.mi), .combine="c") %do% {
    tfi.mi.another[tfi.mi[i,'resnum1'], tfi.mi[i,"resnum2"]]
  }
  
  print('__Input prody MI results:')
  tfi.mip.another <- read.table(file.path(prody.root, stringr::str_c(tfi,"_MIp.txt")), sep='\t', header=F, as.is=T)
  tfi.mi$prody_mip <- foreach(i=1:nrow(tfi.mi), .combine="c") %do% {
    tfi.mip.another[tfi.mi[i,'resnum1'], tfi.mi[i,"resnum2"]]
  }

  
  mi.list[[tfi]] <- tfi.mi
  
}

## calculate the pearson correlation between the MI values from ProDy and MISTIC

p.list_mistic.scaled_prody.mi <- list()
p.list_mistic.scaled_prody.mip <- list()

for(tfi in tf.fams){
  print(tfi)
  
  prody.value <- (mi.list[[tfi]]$prody_value)
  mistic.value <- ( mi.list[[tfi]]$mistic_value)
  mistic.scaled <- foreach(i=mistic.value, .combine='c') %do%{
    (i-min(mistic.value)) / (max(mistic.value)-min(mistic.value))
  }
  
  
  curr.mi.mat <- mi.list[[tfi]]
  curr.mi.mat$mistic_scaled <- mistic.scaled
  
  curr.mi.mat <- curr.mi.mat[curr.mi.mat$mistic_scaled>0, ]
  
  # cor.test(as.numeric( curr.mi.mat$mistic_scaled), as.numeric(curr.mi.mat$prody_value))
  # cor.test(as.numeric( curr.mi.mat$mistic_scaled), as.numeric(curr.mi.mat$prody_mip))
  
  p.list_mistic.scaled_prody.mi[[tfi]] <- ggscatter(curr.mi.mat, x = "prody_value", y = "mistic_scaled", add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE  ) + stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, method="spearman") +  ggtitle(tfi) 
    
    
  p.list_mistic.scaled_prody.mip[[tfi]] <- ggscatter(curr.mi.mat, x = "prody_mip", y = "mistic_scaled", add = "reg.line",  # Add regressin line
                                                     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                                     conf.int = TRUE  ) + stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, method="spearman") +  ggtitle(tfi) 
  
}


ggexport(plotlist = p.list_mistic.scaled_prody.mi, filename = "E:/plots_mistic.scaled.vs.prody.mi_spearman.pdf",
         nrow = 3, ncol = 3)

ggexport(plotlist = p.list_mistic.scaled_prody.mip, filename = "E:/plots_mistic.scaled.vs.prody.mip_spearman.pdf",
         nrow = 3, ncol = 3)


corr.vec <- c(0.73, 0.5, 0.63, .7, .73, .74, .69, .66, .73)
