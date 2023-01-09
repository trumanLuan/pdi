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

