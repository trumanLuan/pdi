## :::::::::::::::::::::::::::::::::::::::::::::::
## required packages
## :::::::::::::::::::::::::::::::::::::::::::::::

require(dplyr)
require(ggplot2)
library(RColorBrewer)
require(reshape2)
require(stringr)

require(MotIV)
require(seqLogo)
require(MotifDb)
require(motifStack)


## :::::::::::::::::::::::::::::::::::::::::::::::
## project dirs and data files.
## :::::::::::::::::::::::::::::::::::::::::::::::

## dirs
proj.dir="E:/projects/PDI"
code.dir <- file.path(proj.dir, "code")
out.dir=file.path(proj.dir, "results")
data.dir = file.path(proj.dir, 'data')
tf2dna.data.dir <- file.path(data.dir, "TF2DNA")
tf2dna.tf.seq.dir <- file.path(data.dir, "TF2DNA/tf_domain_aln")

## data files
f.tf2dna <- file.path(data.dir, "TF2DNA/tf2dna.experiments.unique.csv")
tf2dna=read.csv(f.tf2dna, header=TRUE, as.is=TRUE)

species = c('Homo_sapiens', 'Mus_musculus', 'Drosophila_melanogaster', 'Saccharomyces_cerevisiae', 'Caenorhabditis_elegans' )
tf2dna <- tf2dna %>% subset(TF_Species %in% species)

tf.maj.fam=table(tf2dna[,'DBDs']); 
tf2dna <- tf2dna %>% subset(DBDs %in% names(tf.maj.fam[tf.maj.fam>30]))

load(file.path(data.dir, 'TF2DNA/dna_motifs', 'motif_summary_pwm.RData'))  # load pwm.list, d.stamp, sub.sub.tf2dna


## :::::::::::::::::::::::::::::::::::::::::::::::
## Compare our collected DNA motifs to known databases.
## :::::::::::::::::::::::::::::::::::::::::::::::

compar.motifdb.dt <- read.table(file.path(tf2dna.data.dir, "compar_dna_motifs_to_known_motifdb.txt"), header=T, sep='\t', as.is=T ) %>%
  subset(annot_DB %in% c("HOCOMOCOv11", "jaspar2018", "UniPROBE"))

## how many TFs voerlapped with known ones in jaspar database for each TF family?
motifs.to.annot.uniq <- unique(compar.motifdb.dt[,1:4] )

tf.fam.vec <- names(table(tf2dna$DBDs))
table(motifs.to.annot.uniq$tf_fam )[tf.fam.vec]
a <- table(motifs.to.annot.uniq$tf_fam )[tf.fam.vec]/c(65, 41, 42,123,45,294,198,54,41 )
mean(a)


## compare the similarity of dna motifs between our collected dnaset and those in these databases

compar.motifdb.dt$index <- stringr::str_c(compar.motifdb.dt$species, "+", compar.motifdb.dt$tf_fam, "+", compar.motifdb.dt$tf_name, "+",compar.motifdb.dt$tf_ID)
tf2dna$index <- stringr::str_c(tf2dna$TF_Species, "+", tf2dna$DBDs, "+", tf2dna$TF_Name, "+", tf2dna$DBID)
tf2dna$annotated <- ifelse(tf2dna$index %in% compar.motifdb.dt$index, "Yes", "No")


motif.dist.matrix <- as.matrix(d.stamp)
corr.df.with.annot = NULL
annot.motifs = subset(tf2dna, annotated %in% "Yes")
jaspar.scores <- MotIV::readDBScores(paste(system.file(package="MotIV"),"/extdata/jaspar2010_PCC_SWU.scores",sep=""))
pb = txtProgressBar(min = 0, max = nrow(annot.motifs), initial = 0, style = 3) 

for(i in 1:nrow(annot.motifs)){
  tmp.index = annot.motifs$index[i]
  tmp.motif.id = annot.motifs$Motif_ID[i]
  tmp.annot.motifs = subset(compar.motifdb.dt, index %in% tmp.index)
  
  tmp.dist = sapply(tmp.annot.motifs$annot_motif, function(x) {
    a=MotIV::motifDistances(list(our.pwm=pwm.list[[tmp.motif.id]], annot.pwm=MotifDb[[as.character(x)]] ), DBscores=jaspar.scores, cc="PCC", align="SWU", top=5, go=1, ge=0.5)
    as.matrix(a)[2,1]
  })
  
  tmp.similar = 1-tmp.dist
  tmp.annot.motifs$corr.with.annot = tmp.similar
  corr.df.with.annot = rbind(corr.df.with.annot, tmp.annot.motifs)
  
  setTxtProgressBar(pb,i)
}; close(pb)

corr.df.with.annot$annot_tf_in_DB <- toupper(corr.df.with.annot$annot_tf_in_DB)
corr.df.with.annot$tf_name <- toupper(corr.df.with.annot$tf_name)
table(corr.df.with.annot$tf_name == corr.df.with.annot$annot_tf_in_DB )
is.in <- corr.df.with.annot$tf_name == corr.df.with.annot$annot_tf_in_DB 
test <- corr.df.with.annot; test$is_common <- is.in
write.table(test, file.path(tf2dna.data.dir, "compar_dna_motifs_to_known_motifdb_motifMatch_results.txt"), quote=F,sep='\t', row.names = F, col.names=T)

## hist of similarity score between our collected data and other databases.
ggplot(corr.df.with.annot, aes(x=corr.with.annot)) + 
  geom_histogram() + theme_classic() +
  facet_grid(tf_fam ~ .)

max.corr.df.with.annot = corr.df.with.annot %>% group_by(tf_fam, tf_name, tf_ID, annot_DB) %>% summarise(corr.with.annot=max(corr.with.annot))
ggplot(max.corr.df.with.annot, aes(x=tf_fam, y=corr.with.annot)) + 
  geom_boxplot() + theme_classic() +
  coord_flip() + facet_grid(. ~ annot_DB)

a <- max.corr.df.with.annot %>% group_by(tf_fam) %>% summarise(mean=mean(corr.with.annot))
mean(a$mean)

## visualization of sequence logos of selected dna motifs with seqLogo
seqLogo::seqLogo(MotifDb[["Hsapiens-jaspar2018-NFE2-MA0841.1"]]) ## for bzip
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+bZIP_1+NFE2+ENSG00000123405")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Mmusculus-UniPROBE-Ehf.UP00015"]]) ## for Ets
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+Ets+EHF+ENSG00000135373")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Hsapiens-jaspar2018-FOXO3-MA0157.2"]]) ## for fork head
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+Fork_head+FOXO3+ENSG00000118689")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Hsapiens-JASPAR_2014-USF1-MA0093.2"]]) ## for HLH
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+HLH+USF1+ENSG00000158773")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Hsapiens-jaspar2018-SOX4-MA0867.1"]]) ## for HMG box
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+HMG_box+SOX4+ENSG00000124766")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Hsapiens-JASPAR_CORE-HNF1B-MA0153.1"]]) ## for Hoemmobox
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+Homeobox+HNF1B+ENSG00000108753")$Motif_ID]] )

# seqLogo::seqLogo(MotifDb[["Hsapiens-jaspar2018-YY1-MA0095.1"]]) ## for zf-c2h2
# seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+zf-C2H2+YY1+ENSG00000100811")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Mmusculus-UniPROBE-Zscan4.UP00026"]]) ## for zf-c2h2, ZSCAN4
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Mus_musculus+zf-C2H2+Zscan4+EDL38120.1")$Motif_ID]] )

# seqLogo::seqLogo(MotifDb[["Hsapiens-jaspar2018-VDR-MA0693.2"]]) ## for zf-c4
# seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Homo_sapiens+zf-C4+VDR+ENSG00000111424")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Mmusculus-jaspar2018-Ar-MA0007.3"]]) ## for zf-c4: AR
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Mus_musculus+zf-C4+Ar+ENSMUSG00000046532")$Motif_ID]] )

seqLogo::seqLogo(MotifDb[["Scerevisiae-jaspar2018-UGA3-MA0410.1"]]) ## for zn clus
seqLogo::seqLogo( pwm.list[[subset(annot.motifs, index %in% "Saccharomyces_cerevisiae+Zn_clus+UGA3+YDL170W")$Motif_ID]] )

## :::::::::::::::::::::::::::::::::::::::::::::::
## Motif alignment and merging.
## :::::::::::::::::::::::::::::::::::::::::::::::

tf.fam = (unique(tf2dna$DBDs))
plot.motifs = list()

for(curr.tf.fam in tf.fam){
  
  cat(curr.tf.fam, 'BEGIN\n')
  curr.motifs = readLines(file.path(data.dir, 'TF2DNA/dna_motifs/raw_by_tf_fam', paste0(curr.tf.fam, '_motif_matrix.txt'))) 
  curr.motifs = curr.motifs[grepl('>M', curr.motifs)] 
  curr.motifs = stringr::str_replace_all(curr.motifs, '>', '')
  
  ## motif merge and motifs cluster logos
  sub.pfms = pwm.list[curr.motifs]
  sub.pfms <- sub.pfms[sapply(sub.pfms, is.numeric)]
  sub.pfms <- lapply(names(sub.pfms), function(.ele, sub.pfms){new("pfm",mat=sub.pfms[[.ele]], name=.ele)}, 
                     sub.pfms) 
  
  sub.pfms.aln = motifStack::DNAmotifAlignment(pfms=sub.pfms) # motifs alignment.
  pfmat.aln = lapply(sub.pfms.aln, function(x) x@mat)
  names(pfmat.aln) = sapply(sub.pfms.aln, function(x) x@name)
  pfmat.aln.new = Reduce('+', pfmat.aln)
  pfmat.aln.new = pfmat.aln.new/(length(pfmat.aln))
  motif <- new("pfm", mat=pfmat.aln.new, name=curr.tf.fam)
  #motif <- trimMotif(motif, t=0.2)
  plot.motifs[[curr.tf.fam]] = motif
}

Sys.setenv(R_GSCMD="\"C:/Program Files/gs/gs9.22/bin/gswin64c.exe\"")
pdf(file.path(tf2dna.data.dir, paste('motifs_stack.pdf',sep='_')) )
motifStack::plotMotifLogoStack(plot.motifs[1], ncex=0.8) #
dev.off()

## :::::::::::::::::::::::::::::::::::::::::::::::
## compare the similarity of motifs for each TF family
## :::::::::::::::::::::::::::::::::::::::::::::::

motif.dist.mat = as.matrix(d.stamp)
motif.corr.mat = 1-motif.dist.mat
motif.corr.df = reshape2::melt(motif.corr.mat)
colnames(motif.corr.df) = c("motif1", "motif2", "Corr")

motif.corr.df$tf_fam1 = sapply(motif.corr.df$motif1, function(x) subset(tf2dna, Motif_ID %in% x)$DBDs )
motif.corr.df$tf_fam2 = sapply(motif.corr.df$motif2, function(x) subset(tf2dna, Motif_ID %in% x)$DBDs )
motif.corr.df = subset(motif.corr.df, tf_fam1 == tf_fam2)

annot.motif.vector = subset(tf2dna, annotated %in% "Yes")$Motif_ID
motif.corr.df$annot_motif1 = sapply(motif.corr.df$motif1, function(x) ifelse(x %in% annot.motif.vector, "Yes", "No") )
motif.corr.df$annot_motif2 = sapply(motif.corr.df$motif2, function(x) ifelse(x %in% annot.motif.vector, "Yes", "No") )
motif.corr.df$annot = stringr::str_c(motif.corr.df$annot_motif1, "_", motif.corr.df$annot_motif2)
motif.corr.df$annot[which(motif.corr.df$annot == "No_Yes", arr.ind=T)] = "Yes_No"

## visualization
ggplot(motif.corr.df, aes(x=tf_fam1, y=Corr)) + geom_boxplot() + theme_classic()
# ggplot(motif.corr.df, aes(x=tf_fam1, y=Corr, fill=annot)) + geom_boxplot() + theme_classic()


## visualization, group by different species
motif.corr.df$species_motif1 = sapply(motif.corr.df$motif1, function(x) subset(tf2dna, Motif_ID %in% x)$TF_Species )
motif.corr.df$species_motif2 = sapply(motif.corr.df$motif2, function(x) subset(tf2dna, Motif_ID %in% x)$TF_Species )
motif.corr.df = subset(motif.corr.df, species_motif1 == species_motif2 )

ggplot(motif.corr.df, aes(x=tf_fam1, y=Corr)) + geom_boxplot() + theme_bw()+ facet_grid(species_motif1~.)


## :::::::::::::::::::::::::::::::::::::::::::::::
## whether the variation of between-motif correlation was related to the TF numbers?
## 
## :::::::::::::::::::::::::::::::::::::::::::::::

quantile.range.df = motif.corr.df %>% group_by(species_motif1, tf_fam1) %>% 
  summarise(Q1=quantile(Corr, 0.25, na.rm=T), Q3=quantile(Corr, 0.75, na.rm=T))
quantile.range.df$IQR = quantile.range.df$Q3 - quantile.range.df$Q1

motif.count.by.species = tf2dna %>% group_by(TF_Species, DBDs) %>% summarise(n.motif=n())

all(quantile.range.df$species_motif1 == motif.count.by.species$TF_Species)
all(quantile.range.df$tf_fam1 == motif.count.by.species$DBDs)
quantile.range.df$n.motif = motif.count.by.species$n.motif

## visualization
ggplot(quantile.range.df, aes(x=IQR, y=log2(n.motif),  color=tf_fam1)) +  
  geom_point(aes(col=species_motif1))+
  theme_bw() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+facet_wrap(.~tf_fam1)

## statistical significance testing.
tf.vector = unique(quantile.range.df$tf_fam1)
cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[1])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[1])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[2])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[2])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[3])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[3])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[4])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[4])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[5])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[5])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[6])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[6])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[7])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[7])$n.motif, method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[8])$IQR, 
         log2(subset(quantile.range.df, tf_fam1 %in% tf.vector[8])$n.motif), method="spearman")

cor.test(subset(quantile.range.df, tf_fam1 %in% tf.vector[9])$IQR, 
         subset(quantile.range.df, tf_fam1 %in% tf.vector[9])$n.motif, method="spearman")


