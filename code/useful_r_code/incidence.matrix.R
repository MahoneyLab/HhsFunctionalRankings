#This function creates an incidence matrix from a two-column
#matrix for creating bi-partite graphs


incidence.matrix <- function(mat, verbose = FALSE){
	
	var1 <- unique(mat[,1])
	var2 <- unique(mat[,2])
	inc.mat <- matrix(0, nrow = length(var1), ncol = length(var2))
	rownames(inc.mat) <- var1
	colnames(inc.mat) <- var2
	
	if(length(var1) <= length(var2)){
		for(i in 1:length(var1)){
			if(verbose){report.progress(i, length(var1))}
			var1.locale <- which(mat[,1] == var1[i])
			var2.links <- mat[var1.locale,2]
			var2.locale <- match(var2.links, var2)
			test <- which(var2 %in% var2.links)
			inc.mat[i,var2.locale] <- 1
		}
	}
	
	if(length(var2) < length(var1)){
		for(i in 1:length(var2)){
			if(verbose){report.progress(i, length(var2))}
			var2.locale <- which(mat[,2] == var2[i])
			var.links <- mat[var2.locale,1]
			var1.locale <- match(var.links, var1)
			inc.mat[var1.locale,i] <- 1
		}
	}
	
	cat("\n")
	return(inc.mat)
}