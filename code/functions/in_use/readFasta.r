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

