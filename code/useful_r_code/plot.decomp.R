#This function mean centers a matrix,
#takes the svd and then plots the specified
#principal components against each other.
#cols is a vector of color values for the points
#color.by is a numeric matrix of factors whose
#relationship to the principal components can 
#be examined. 
#color.by supercedes cols if specified.

plot.decomp <- function(mat, pc = 2, cols = rep("black", nrow(mat)), color.by = NULL, pch = NULL, main = NULL){

	cent <- apply(mat, 2, function(x) mean(x) - x)
	decomp <- svd(cent)
	pc.mat <- decomp$u
	colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

	if(is.null(pch)){pch = 16}

	if(is.null(color.by)){
		num.plots = 1
		plot.label <- main
	}else{
		num.plots <- ncol(color.by)
		if(is.null(colnames(color.by))){
			colnames(color.by) <- paste(main, "Factor", 1:ncol(color.by))
		}
		plot.label <- colnames(color.by)
	}
		
	for(i in 1:num.plots){
		if(!is.null(color.by)){
			cols <- colors.from.values(color.by[,i])
		}
		if(pc == 2){
			plot(pc.mat[,1], pc.mat[,2], xlab = "PC1", ylab = "PC2", 
			main = plot.label[i], col = cols, pch = pch)
		}else{
	  		pairs(pc.mat[,1:pc], col = cols, pch = pch, 
			main = plot.label[i])
			}
	}

	invisible(pc.mat)

}
