readFastaAsList<- function(file=file.name) {

data <- readLines(file) # reads file line-wise into vector
if(length(data) == 0) cat("NOTE: Your input file is empty:(\n")

# records.row is a matrix having the start and stop index for each record
records.row.id=grep(">",data)
#name.row=sub(">","",data[records.row.id])


if(length(name.row)==1){
	data.row.id=data.frame(start=2,end=length(data))
} else {
	data.row.id=data.frame(
		start=records.row.id[1:length(records.row.id)]+1,
		end=c(records.row.id[2:length(records.row.id)]-1,length(data)))
}

fasta.lst=list()
for(i in records.row.id){
	curr.name = sub(">","",data[i])
	curr.seq = data[data.row.id[which(records.row.id==i), 1]: data.row.id[which(records.row.id==i), 2]]
	fasta.lst[[curr.name]] = curr.seq
}

return(fasta.lst)
}

