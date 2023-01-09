
##############################
#~ pairs function
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{ usr <- par("usr"); on.exit(par(usr));par(usr = c(0, 1, 0, 1));r <- cor(x, y, use='complete.obs');
	txt <- format(c(r, 0.123456789), digits=digits)[1]; txt <- paste(prefix, txt, sep="");
  if(missing(cex.cor)) cex <- 1;
  test <- cor.test(x,y, use='complete.obs') 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  text(0.5, 0.5, txt, cex = cex);text(.8, .8, Signif, cex=cex, col=2)}

panel.hist <- function(x, ...)
{ usr <- par("usr"); on.exit(par(usr)); par(usr = c(usr[1:2], 0, 1.5) );h <- hist(x, plot = FALSE);
  breaks <- h$breaks; nB <- length(breaks);y <- h$counts; y <- y/max(y);rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)}

panel.smooth <- function (x, y) {points(x, y, pch='.');abline(lm(y~x), col="red")}

##############################
## read fasta file
readFasta<- function(file=file.name) {

data <- readLines(file) # reads file line-wise into vector

# records.row is a matrix having the start and stop index for each record
records.row.id=grep(">",data)
name.row=sub(">","",data[records.row.id])


if(length(name.row)==1){
	data.row.id=data.frame(start=2,end=length(data))
} else {
	data.row.id=data.frame(
		start=records.row.id[1:length(records.row.id)]+1,
		end=c(records.row.id[2:length(records.row.id)]-1,length(data)))
}

data.row=apply(data.row.id,1,function(x)
	{y=data[x[1]:x[2]]
	return(paste(y,collapse=""))})

return(cbind(name.row,data.row))
}

##############################

## transform to fasta format
ToFasta<-function(seq.mat,names=c()){
# -> seq.mat is sequence vector.
# -> names is a character vector, corresponding to the elements in seq.mat.

output=c()

for (i in 1:length(seq(1:length(seq.mat)))){
if(length(names)>0){temp=rbind(paste(">",names[i],sep=""),seq.mat[i])
output=rbind(output,temp)} else
{temp=rbind(paste(">seq",i,sep=""),seq.mat[i])
output=rbind(output,temp)} 

}

return(output)
}
##############################
## read tsv file
read_file = function(file.dir,h=F){
# h: logical value indicating header is or F.
	read.table(file.dir, sep='\t', as.is=T, header=h)
}

##############################
## between-residues distance within the same domain.
betweenResidueDist = function(pdb, tf.chain){
	# -> pdb: object of read.pdb().
	# -> tf.chain: TF domain chain ID in atoms matrix of pdb.
	require(Rpdb)
	sec1=pdb$atoms$chainid == tf.chain
	sec2=pdb$atoms$chainid == tf.chain
	d = distances(pdb, sec1, sec2)
	d = norm(d, type='xyz')
	rownames(d) = pdb$atoms[sec1,'resid']
	colnames(d) = pdb$atoms[sec2,'resid']
	d.new = matrix(0,nrow=length(unique(rownames(d))),ncol=length(unique(colnames(d))))
	rownames(d.new)=unique(rownames(d))
	colnames(d.new)=unique(colnames(d))
	for(m in unique(rownames(d))){
		#m.residue = unique(rownames(d))[m]
		for(n in unique(colnames(d))){
			tmp.d = d[which(rownames(d)==m),which(rownames(d)==n)]
			tmp.d = min(tmp.d)
			d.new[m,n]=tmp.d
	}
}
d.new
}

##############################
## heatmap with data frame or matrix.
my_heatmap = function(gene.exp,main, outputfile){
	library(gplots)
	pdf(outputfile)
	my_palette = colorRampPalette(c("green4","white","deeppink3"))(100)
	# gene clustering
	scale.data = scale(t(gene.exp))
	row.dist = dist(t(scale.data), method='manhattan')
	row.cluster = hclust(row.dist, method='ward')

	heatmap.2(gene.exp, main=main,
          notecol='black', density.info='none', trace='none',
          margins=c(10,9), col=my_palette, dendrogram='row',
          Colv='NA', Rowv=as.dendrogram(row.cluster),
          scale='none', #RowSideColors = row.color,
      cexRow=0.5)
	dev.off()
}





