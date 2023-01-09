## :::::::::::::::::::::::::::::::::::::::::::::::
## required packages
## :::::::::::::::::::::::::::::::::::::::::::::::

rm(list=ls())
require(dplyr)
require(ggplot2)
require(ggpubr)
require(stringr)
require(Rpdb)
require(foreach)


## :::::::::::::::::::::::::::::::::::::::::::::::
## project dirs and data files.
## :::::::::::::::::::::::::::::::::::::::::::::::

## dirs
proj.dir="E:/projects/PDI"
code.dir <- file.path(proj.dir, "code")

out.dir=file.path(proj.dir, "results")
data.dir = file.path(proj.dir, 'data')
pdb.data.dir <- file.path(data.dir, 'pfam_pdb')

## source functions
source( file.path(code.dir, 'functions/in_use/tf2dnaDist.R') )
source( file.path(code.dir, "functions/in_use/readFasta.r") )

## dataset meta and data files
pdb.meta = read.table(file.path(pdb.data.dir,'pfam_pdb_all.txt'),header=T,sep='\t',as.is=T)
pdb.meta = subset(pdb.meta, ! TF_family %in% c('zf-C2H2'))

pdbs.in <- dir(file.path(pdb.data.dir, "pdb"))
pdbs.in <- stringr::str_replace_all(pdbs.in, "\\.pdb", "")
pdb.meta <- subset(pdb.meta, PDB_ID %in% pdbs.in)
pdb.meta$resolution <- rep(NA, nrow(pdb.meta))

length(table(pdb.meta$PDB_ID)) 

## add resolution information into pdb metatable
pdb.fs <- dir(file.path(pdb.data.dir, "pdb"))
for(fi in pdb.fs){
  cat(fi,'\n')
  x.pdb = Rpdb::read.pdb(file.path(pdb.data.dir, "pdb", fi))
  x.pdb.remark = x.pdb$remark
  x.pdb.resolution = grep("REMARK[[:space:]]+2 RESOLUTION", x.pdb.remark, value=T)
  x.pdb.resolution2 = stringr::str_extract(x.pdb.resolution, "\\d+\\.\\d+ ANGSTROMS")
  x.pdb.resolution2 = as.numeric(stringr::str_extract(x.pdb.resolution2, "\\d+\\.\\d+"))
  
  fi.label <- stringr::str_replace(fi, "\\.pdb", "")
  pdb.meta$resolution[ which(pdb.meta$PDB_ID == fi.label, arr.ind=T)  ] <- x.pdb.resolution2
}

pdb.meta.2a <- pdb.meta %>% subset(resolution <=2.5)
length(unique(pdb.meta.2a$PDB_ID))

## SPEER resulting data, for TF subclass determining sites analysis.
load(file.path(data.dir, "intermediate_dt/speer.score.c.rdata"))
speer.score <- speer.score.c; rm(speer.score.c)

## Amino acid alphabet, abbra. of amino acids
aa.alphabet <- read.table(file.path(data.dir,'intermediate_dt/amino_acid_alphabet.txt'), header=T,sep='\t',as.is=T)

## highly coevolving residue paires
hcr.df <- read.table(file.path(data.dir, "intermediate_dt/highly_coevolving_reside_pairs.txt"), header=T,sep='\t',as.is=T )


## :::::::::::::::::::::::::::::::::::::::::::::::
## calculate spatial distance between highly coevolving residue pairs (HCRs)
## SDS were defined with SPEER
## HCRs were defined with the ensemble results from four methods
## :::::::::::::::::::::::::::::::::::::::::::::::

## load TF domain between-residue 3D distance
spatial.resi.dist.f <- list.files(file.path(pdb.data.dir, "pdb_aln"), recursive = F, full.names = F)
spatial.resi.dist.f <-  grep("between.residues.within.domain.distance", spatial.resi.dist.f, value = T)

spatial.resi.dist <- NULL
for(fi in spatial.resi.dist.f){
  message(fi)
  tmp <- read.table(file.path(pdb.data.dir, "pdb_aln", fi), header=T, sep='\t', as.is=T)
  colnames(tmp) <- c( "position1", "position2", "dist_A" )
  tmp$tf_fam <- rep(stringr::str_split(fi, "\\.")[[1]][1], nrow(tmp) )
  spatial.resi.dist <- rbind(spatial.resi.dist, tmp)
}

spatial.resi.dist <- spatial.resi.dist[!spatial.resi.dist$position1 == spatial.resi.dist$position2,]
spatial.resi.dist$resi_pair <- stringr::str_c(spatial.resi.dist$position1, "_", spatial.resi.dist$position2 )

## load TF-DNA between residue 3D distance
spatial.tf2dna.dist.f <- list.files(file.path(pdb.data.dir, "pdb_aln"), recursive = F, full.names = F)
spatial.tf2dna.dist.f <- grep("tf-dna.distance", spatial.tf2dna.dist.f, value=T) 

spatial.tf2dna.dist <- NULL
for(fi in spatial.tf2dna.dist.f){
  message(fi)
  tmp <- read.table(file.path(pdb.data.dir, "pdb_aln", fi), header=F, sep='\t', as.is=T)
  colnames(tmp) <- c( "position", "dist_A" )
  tmp$tf_fam <- rep(stringr::str_split(fi, "\\.")[[1]][1], nrow(tmp) )
  spatial.tf2dna.dist <- rbind(spatial.tf2dna.dist, tmp)
}


## :::::::::::::::::::::::::::::::::::::
## compare TF domain between-residue distance between highly co-evolving residue pairs (HCRs) and non-HCRs
## related to figure_S7.

hcr.df$resi_pair_index <- stringr::str_c(hcr.df$tf_fam, "_", hcr.df$pairs)
spatial.resi.dist$resi_pair_index <- stringr::str_c(spatial.resi.dist$tf_fam, "_", spatial.resi.dist$resi_pair)
spatial.resi.dist$HCR <- ifelse( spatial.resi.dist$resi_pair_index %in% hcr.df$resi_pair_index, "HCR", "nonHCR" )

ggplot2::ggplot(spatial.resi.dist, aes(x=tf_fam, y=dist_A, fill=HCR)) + geom_boxplot() + theme_classic()

## statistical testing
t.test( subset(spatial.resi.dist, tf_fam %in% "bZIP_1" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "bZIP_1" & HCR %in% "nonHCR" )$dist_A )

t.test( subset(spatial.resi.dist, tf_fam %in% "Ets" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "Ets" & HCR %in% "nonHCR" )$dist_A )

t.test( subset(spatial.resi.dist, tf_fam %in% "HLH" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "HLH" & HCR %in% "nonHCR" )$dist_A )

t.test( subset(spatial.resi.dist, tf_fam %in% "HMG_box" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "HMG_box" & HCR %in% "nonHCR" )$dist_A )

t.test( subset(spatial.resi.dist, tf_fam %in% "Homeobox" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "Homeobox" & HCR %in% "nonHCR" )$dist_A )

t.test( subset(spatial.resi.dist, tf_fam %in% "zf-C4" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "zf-C4" & HCR %in% "nonHCR" )$dist_A )

t.test( subset(spatial.resi.dist, tf_fam %in% "Zn_clus" & HCR %in% "HCR" )$dist_A, 
        subset(spatial.resi.dist, tf_fam %in% "Zn_clus" & HCR %in% "nonHCR" )$dist_A )


## count the HCR pairs with length >10A
count.hcr.in.tf.fam <- spatial.resi.dist %>% subset( HCR == "HCR" ) %>% group_by(tf_fam) %>% summarize(count=n())
count.dist.hcr.in.tf.fam <- spatial.resi.dist %>% subset(HCR=="HCR") %>% subset(dist_A >10 ) %>% group_by(tf_fam) %>% summarize(count=n())
mean(count.dist.hcr.in.tf.fam$count / count.hcr.in.tf.fam$count)
count.dist.hcr.in.tf.fam

## :::::::::::::::::::::::::::::::::::::
## In HCRs, compare the between-residues distance between SDS, SDS-nonSDS and intra-nonSDS

dist.between.hcr <- subset(spatial.resi.dist, HCR %in% "HCR")

dist.between.hcr.new <- NULL
for(tfi in unique(dist.between.hcr$tf_fam)){
  message(tfi)
  tmp.speer <- subset(speer.score, tf_fam %in% tfi)
  tmp.dist.hcr <- subset(dist.between.hcr, tf_fam %in% tfi)
  
  sds.site <- as.numeric(subset(tmp.speer, pval < 0.1)$column) + 1 # TSDS sites
  
  site.group = c()
  for(i in 1:nrow(tmp.dist.hcr)){
    if(length(which(tmp.dist.hcr[i,1:2] %in% sds.site)) == 0) site.group = c(site.group, 'extra')
    if(length(which(tmp.dist.hcr[i,1:2] %in% sds.site)) == 1) site.group = c(site.group, 'between')
    if(length(which(tmp.dist.hcr[i,1:2] %in% sds.site)) == 2) site.group = c(site.group, 'intra')
  }
  tmp.dist.hcr$site_group = site.group
  
  dist.between.hcr.new <- rbind(dist.between.hcr.new, tmp.dist.hcr)
  
}

ggplot(dist.between.hcr.new, aes(x=tf_fam, y=dist_A, fill=site_group)) + geom_boxplot() + theme_classic() 

## statistical testing
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "bZIP_1" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "bZIP_1" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "bZIP_1" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "bZIP_1" & site_group %in% "intra" )$dist_A )

wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "Ets" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "Ets" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "Ets" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "Ets" & site_group %in% "intra" )$dist_A )

wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "HLH" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "HLH" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "HLH" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "HLH" & site_group %in% "intra" )$dist_A )

wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "HMG_box" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "HMG_box" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "HMG_box" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "HMG_box" & site_group %in% "intra" )$dist_A )

wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "Homeobox" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "Homeobox" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "Homeobox" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "Homeobox" & site_group %in% "intra" )$dist_A )

wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "zf-C4" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "zf-C4" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "zf-C4" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "zf-C4" & site_group %in% "intra" )$dist_A )

wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "Zn_clus" & site_group %in% "between" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "Zn_clus" & site_group %in% "intra" )$dist_A )
wilcox.test( subset(dist.between.hcr.new, tf_fam %in% "Zn_clus" & site_group %in% "extra" )$dist_A, 
             subset(dist.between.hcr.new, tf_fam %in% "Zn_clus" & site_group %in% "intra" )$dist_A )


## :::::::::::::::::::::::::::::::::::::
## Compare the distance to DNA interface of co-evolving residues: two residues in HCRs

dist.tf2dna.hcr <- NULL
for(tfi in unique(spatial.tf2dna.dist$tf_fam)){
  message(tfi)
  ## the position closer to DNA interface defined as position1
  ## the position more distant from DNA interface defined as position2
  
  tmp.dist.in.hcr <- subset(dist.between.hcr.new, tf_fam %in% tfi)
  
  for(i in 1:nrow(tmp.dist.in.hcr)){
    tmp.p1.tf2dna <- subset(spatial.tf2dna.dist, tf_fam %in% tfi & position %in% tmp.dist.in.hcr[i,"position1"])$dist_A
    tmp.p2.tf2dna <- subset(spatial.tf2dna.dist, tf_fam %in% tfi & position %in% tmp.dist.in.hcr[i,"position2"])$dist_A
    
    if(any(is.na(tmp.p1.tf2dna) | is.na(tmp.p2.tf2dna))) next;
    if(tmp.p1.tf2dna < tmp.p2.tf2dna){
      to.p1 <- tmp.dist.in.hcr[i,"position1"]
      to.p2 <- tmp.dist.in.hcr[i,"position2"]
      tf2dna.p1 <- tmp.p1.tf2dna
      tf2dna.p2 <- tmp.p2.tf2dna
    }else{
      to.p1 <- tmp.dist.in.hcr[i,"position2"]
      to.p2 <- tmp.dist.in.hcr[i,"position1"]
      tf2dna.p1 <- tmp.p2.tf2dna
      tf2dna.p2 <- tmp.p1.tf2dna
    }
    
    tmp.output.df <- data.frame(tf_fam=tfi, resi_pairs=tmp.dist.in.hcr[i, "resi_pair"], tsds_site_group= tmp.dist.in.hcr[i,"site_group"], 
                                position1=to.p1, position2=to.p2, p1.dist2dna=tf2dna.p1, p2.dist2dna=tf2dna.p2, stringsAsFactors = F  )
    dist.tf2dna.hcr <- rbind(dist.tf2dna.hcr, tmp.output.df)
  }
}

## plot
ggplot(dist.tf2dna.hcr, aes(x=p1.dist2dna, y=p2.dist2dna, col=tsds_site_group)) + geom_point(size=2, shape=18) +
  facet_wrap(.~tf_fam) + scale_color_brewer(palette="Set1") +
  # geom_smooth(method=lm,  linetype="dashed", color="darkred", fill="blue") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=1) +
  xlim(0,30) + ylim(0, 30) + theme_bw()

# write.table(dist.tf2dna.hcr, file=file.path(data.dir, "dist_tf2dna_of_HCR.txt"), quote=F,sep='\t', row.names = F, col.names = T)

## add TF-domain between-residue distance information into dist.tf2dna.hcr.
inter.resi.dist <- NULL
for(i in 1:nrow(dist.tf2dna.hcr)){
  tmp.dist <- subset(dist.between.hcr, tf_fam == dist.tf2dna.hcr$tf_fam[i] & resi_pair == dist.tf2dna.hcr$resi_pairs[i])$dist_A
  inter.resi.dist <- c(inter.resi.dist, tmp.dist)
}
dist.tf2dna.hcr$inter.resi.dist <- inter.resi.dist

## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for Homeobox
a <- subset(dist.tf2dna.hcr, tf_fam == "Homeobox")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()


## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for HLH
a <- subset(dist.tf2dna.hcr, tf_fam == "HLH")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()



## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for Ets
a <- subset(dist.tf2dna.hcr, tf_fam == "Ets")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()


## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for bZIP_1
a <- subset(dist.tf2dna.hcr, tf_fam == "bZIP_1")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for HMG_box
a <- subset(dist.tf2dna.hcr, tf_fam == "HMG_box")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()


## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for zf-C4
a <- subset(dist.tf2dna.hcr,tf_fam == "zf-C4")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA
table(a$is.label) # 32.6%

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()


## ->
## scatter plot for comparing the dist to DNA interface of each residue in CRP
## for Zn_clus

a <- subset(dist.tf2dna.hcr, tf_fam == "Zn_clus")
a$is.label <- ifelse(a$p1.dist2dna>10 | a$p2.dist2dna>10, "YES", "NO")
a$pairs.label <- a$resi_pairs
a$pairs.label[which(a$is.label == "NO", arr.ind = T)] <- NA
table(a$is.label) # 33.3%

library(ggpubr)
ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07"), label = "pairs.label", repel = TRUE
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()

ggscatter(
  a, x = "p1.dist2dna", y = "p2.dist2dna",
  color = "tsds_site_group", size = "inter.resi.dist", alpha = 0.6, 
  palette = c("#00AFBB", "#E7B800", "#FC4E07")
) + geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred", size=0.5)+
  geom_hline(yintercept = 10, linetype="dashed", color="darkred", size=0.5) +
  geom_vline(xintercept = 10, linetype="dashed", color="darkred", size=0.5) + theme_bw()


## the average percentage of residues with a distance to DNA interface > 10A
p.residue.dist.dna=c(0.281, 0.451, 0.536, 0.675, 0.186, 0.326, .333)

## :::::::::::::::::::::::::::::::::::::::::::::::
## for revised manuscript
## aim: compare the spatial distance between tf domain residues and between tf-dna residues
## 
## :::::::::::::::::::::::::::::::::::::::::::::::

## files
dist.tf2tf.4a.fs <- grep( "between.residues.within.domain.distance", dir(file.path(pdb.data.dir, "pdb_aln")), value=T )
dist.tf2dna.4a.fs <- grep("tf-dna.distance", dir(file.path(pdb.data.dir, "pdb_aln")), value=T)

dist.tf2tf.2.5a.fs <- grep( "between.residues.within.domain.distance", dir(file.path(pdb.data.dir, "pdb_2a_reso")), value=T )
dist.tf2dna.2.5a.fs <-  grep("tf-dna.distance", dir(file.path(pdb.data.dir, "pdb_2a_reso")), value=T)

## load data content
dist.tf2tf.4a <- foreach(fi = dist.tf2tf.4a.fs, .combine="rbind") %do% {
  message(fi)
  tfi <- stringr::str_split(fi, "\\.")[[1]][1]
  
  test.dist.tf2tf <- read.table(file.path(pdb.data.dir, "pdb_aln", fi ), header=T,sep='\t', as.is=T)
  colnames(test.dist.tf2tf) <- c("p1", 'p2', 'dist')
  test.dist.tf2tf$resi_pair <- stringr::str_c(test.dist.tf2tf$p1, '_', test.dist.tf2tf$p2)
  test.dist.tf2tf$tf_fam <- rep(tfi, nrow(test.dist.tf2tf))
  test.dist.tf2tf
}

dist.tf2tf.2.5a <- foreach(fi = dist.tf2tf.2.5a.fs, .combine="rbind") %do% {
  message(fi)
  tfi <- stringr::str_split(fi, "\\.")[[1]][1]
  
  test.dist.tf2tf <- read.table(file.path(pdb.data.dir, "pdb_2a_reso", fi ), header=T,sep='\t', as.is=T)
  colnames(test.dist.tf2tf) <- c("p1", 'p2', 'dist')
  test.dist.tf2tf$resi_pair <- stringr::str_c(test.dist.tf2tf$p1, '_', test.dist.tf2tf$p2)
  test.dist.tf2tf$tf_fam <- rep(tfi, nrow(test.dist.tf2tf))
  test.dist.tf2tf
}

dist.tf2dna.4a <- foreach(fi = dist.tf2dna.4a.fs, .combine="rbind") %do% {
  message(fi)
  tfi <- stringr::str_split(fi, "\\.")[[1]][1]
  
  test.dist.tf2dna <- read.table(file.path(pdb.data.dir, "pdb_aln", fi ), header=F,sep='\t', as.is=T)
  colnames(test.dist.tf2dna) <- c("position", 'dist')
  test.dist.tf2dna$tf_fam <- rep(tfi, nrow(test.dist.tf2dna))
  test.dist.tf2dna
}

dist.tf2dna.2.5a <- foreach(fi = dist.tf2dna.2.5a.fs, .combine="rbind") %do% {
  message(fi)
  tfi <- stringr::str_split(fi, "\\.")[[1]][1]
  
  test.dist.tf2dna <- read.table(file.path(pdb.data.dir, "pdb_2a_reso", fi ), header=F,sep='\t', as.is=T)
  colnames(test.dist.tf2dna) <- c("position", 'dist')
  test.dist.tf2dna$tf_fam <- rep(tfi, nrow(test.dist.tf2dna))
  test.dist.tf2dna
}

## compare the similarity between two versions of spatial distance measurements.
dist.tf2tf.4a$resolution <- rep("4a", nrow(dist.tf2tf.4a))
dist.tf2tf.2.5a$resolution <- rep("2.5a", nrow(dist.tf2tf.2.5a))
dist.tf2tf <- dist.tf2tf.4a; 
if(all(dist.tf2tf.4a$position == dist.tf2tf.2.5a$position)){
  dist.tf2tf$dist_4a <- dist.tf2tf$dist;
  dist.tf2tf$dist_2.5a <- dist.tf2tf.2.5a$dist
}

dist.tf2dna.4a$resolution <- rep("4a", nrow(dist.tf2dna.4a))
dist.tf2dna.2.5a$resolution <- rep("2.5a", nrow(dist.tf2dna.2.5a))
dist.tf2dna <- dist.tf2dna.4a;
if(all(dist.tf2dna.4a$position == dist.tf2dna.2.5a$position)){
  dist.tf2dna$dist_4a <- dist.tf2dna.4a$dist
  dist.tf2dna$dist_2.5a <- dist.tf2dna.2.5a$dist
}


ggplot(dist.tf2tf, aes(x = dist_4a, y=dist_2.5a)) + geom_point(size=2, alpha=0.6) +
  geom_smooth(method=lm, se=FALSE) +
  facet_wrap(~tf_fam)

cor.test( x=subset(dist.tf2tf, tf_fam == "bZIP_1")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'bZIP_1')$dist_2.5a  )
cor.test( x=subset(dist.tf2tf, tf_fam == "Ets")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'Ets')$dist_2.5a  )
cor.test( x=subset(dist.tf2tf, tf_fam == "HLH")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'HLH')$dist_2.5a  )
cor.test( x=subset(dist.tf2tf, tf_fam == "HMG_box")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'HMG_box')$dist_2.5a  )
cor.test( x=subset(dist.tf2tf, tf_fam == "Homeobox")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'Homeobox')$dist_2.5a  )
cor.test( x=subset(dist.tf2tf, tf_fam == "zf-C4")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'zf-C4')$dist_2.5a  )
cor.test( x=subset(dist.tf2tf, tf_fam == "Zn_clus")$dist_4a, y=subset(dist.tf2tf, tf_fam == 'Zn_clus')$dist_2.5a  )


ggplot(dist.tf2dna, aes(x = dist_4a, y=dist_2.5a)) + geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=FALSE) +
  facet_wrap(~tf_fam)

cor.test( x=subset(dist.tf2dna, tf_fam == "bZIP_1")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'bZIP_1')$dist_2.5a  )
cor.test( x=subset(dist.tf2dna, tf_fam == "Ets")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'Ets')$dist_2.5a  )
cor.test( x=subset(dist.tf2dna, tf_fam == "HLH")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'HLH')$dist_2.5a  )
cor.test( x=subset(dist.tf2dna, tf_fam == "HMG_box")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'HMG_box')$dist_2.5a  )
cor.test( x=subset(dist.tf2dna, tf_fam == "Homeobox")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'Homeobox')$dist_2.5a  )
cor.test( x=subset(dist.tf2dna, tf_fam == "zf-C4")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'zf-C4')$dist_2.5a  )
cor.test( x=subset(dist.tf2dna, tf_fam == "Zn_clus")$dist_4a, y=subset(dist.tf2dna, tf_fam == 'Zn_clus')$dist_2.5a  )


## :::::::::::::::::::::::::::::::::::::::::::::::
## computational mutation analyses using FoldX.
## :::::::::::::::::::::::::::::::::::::::::::::::



