#This function ranks elements of a matrix based on two
#dimensions. The rank of each element is based on how many 
#elements are contained in the square determined by the 
#element's x and y coordinates. If many elements are in this
#square, the element is highly ranked. If only a few elements
#are in this square, the element is ranked low.
#plot(mat[,1], mat[,2], type = "n");text(mat[,1], mat[,2], labels = rownames(mat))

rank.2D <- function(mat, return.ordered = FALSE){

    get.rank <- function(idx){
        if(any(is.na(mat[idx,]))){return(NA)}

        above.x <- which(mat[,1] >= mat[idx,1])
        above.y <- which(mat[,2] >= mat[idx,2])
        above.both <- length(intersect(above.x, above.y))

        #col <- rep("black", nrow(mat))
        #col[intersect(above.x, above.y)] <- "red"
        #plot(mat[,1], mat[,2], col = col, pch = 16)
        #abline(v = mat[idx,1])
        #abline(h = mat[idx,2])

        return(above.both)
    }

    rank.mat <- matrix(sapply(1:nrow(mat), get.rank), ncol = 1)
    rownames(rank.mat) <- rownames(mat)
    colnames(rank.mat) <- "rank"
    ordered.rank <- rank.mat[order(rank.mat[,1], decreasing = FALSE),,drop=FALSE]
 
    if(return.ordered){
        return(ordered.rank)
    }else{
        return(rank.mat)
    }
}