#This script subsets a matrix using an expression. 
#The expression is in the following form:
#column.name == value, or column.name != value
#or calumn.name < value, etc.
#where colum.name is an unambiguous string
#identifying a particular column
#The script returns the subsetted cross.

subset.expr <- function(mat, expr){
	
	sub.mat <- mat
	
	expr.pieces <- strsplit(expr, "\\ ")
	if(length(expr.pieces[[1]]) != 3){
		stop("Expression must be in the format 'colname comparison value'")
		}
	
	col.locale <- grep(as.character(expr.pieces[[1]][1]), colnames(sub.mat))
	if(length(col.locale) == 0){
		stop("I can't find the column name: ", expr.pieces[[1]][1])
		}
		
	if(length(col.locale) > 1){
		stop("There is more than one column that matches the string: ", expr.pieces[[1]][1])
		}

	#if the final element in the expression is a number
	if(!is.na(as.numeric(expr.pieces[[1]][3]))){
		cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,col.locale]), as.numeric(expr.pieces[[1]][3]))
		}else{ #otherwise
			cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,col.locale]), expr.pieces[[1]][3])
			}
			
	vals.locale <- which(eval(cl))	
	
	sub.mat <- mat[vals.locale,]
	return(sub.mat)
	
}