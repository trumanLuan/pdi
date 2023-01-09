list2matrix = function(list, byrow=F, bycol=F){
# transform your list to a matrix. firstly, you have a list that contain a series of matrix in each label, these matrix have a same number of rows or columns.
# -> byrow: matrix have same rows. call cbind
# -> bycol: matrix have same columns. call rbind.

matrix = NULL
if(byrow & bycol) cat("ERROR: You need to combine sub-matrix by row, or by col exclusively.\n")

if(byrow){
	for(i in 1:length(list)){
		matrix = cbind(matrix, list[[i]])
	}
	return(matrix)
}

if(bycol){
	for(i in 1:length(list)){
		matrix = rbind(matrix, list[[i]])
	}
	return(matrix)
}

}
