#This function uses a Fisher's Exact test to test 
#overlap between two vectors.

test.overlap <- function(V1, V2, all.V){

    set1.size <- length(V1)
    set2.size <- length(V2)
    overlap.size <- length(intersect(V1, V2))
    set1.not2 <- length(setdiff(V1, V2))
    set2.not1 <- length(setdiff(V2, V1))
    either <- union(V1, V2)
    neither <- length(setdiff(all.V, either))
    

    test.mat <- matrix(c(neither, set2.not1, set1.not2, overlap.size), nrow = 2, byrow = TRUE)
    colnames(test.mat) <- rownames(test.mat) <- c("no", "yes")

    f.p <- fisher.test(test.mat, alternative = "two.sided")$p.value
    results <- list("contingency.table" = test.mat, "p" = f.p)
    return(results)
}