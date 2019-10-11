#This function clusters a network iteratively with fast greedy
#clustering. It uses module size as a stopping criterion.
#To suppress iterative clustering, set max.mod.size to NULL
#If a min.mod.size is set, no modules smaller than that size
#will be produced. If modules smaller than the miminum size 
#are produced they are merged together into one "garbage module" 
#as they are in WGCNA .
#merging cluster can re-generate clusters larger
#than the max, which should be taken into account when setting
#both the min and max.

iter.cluster2 <- function(net, max.mod.size = 400, min.mod.size = 100, 
plot.cluster.size = FALSE, merge.small.clusters = FALSE, sep = "-"){
	
	if(is.null(V(net)$name)){V(net)$name <- 1:vcount(net)}

	clust <- cluster_fast_greedy(net)
	mods <- clust$membership
	
	if(is.null(max.mod.size)){
		return(mods)
		}

	#==========================================
	# internal functions
	#==========================================
	make.clust.tree <- function(clust.obj){
		clust.tree <- vector(mode = "list", length = length(clust.obj))
		for(i in 1:length(clust.obj)){
			clust.tree[[i]] <- clust.obj[[i]]
			}
		return(clust.tree)
		}
				
	get.mod.size <- function(clust.branch){
		if(class(clust.branch) != "list"){return(length(clust.branch))}
		if(class(clust.branch) == "list"){lapply(clust.branch, get.mod.size)}	
		}	
	
	subcluster <- function(clust.branch){

		if(class(clust.branch) != "list"){
			mod.size <- length(clust.branch)
			if(mod.size > max.mod.size){
				subnet <- induced_subgraph(net, vids = clust.branch)
				V(subnet)$name <- clust.branch
				sub.clust <- cluster_fast_greedy(subnet)
				sub.branch <- make.clust.tree(sub.clust)
				return(sub.branch)
			}else{
				return(clust.branch)
			}
		}
			if(class(clust.branch) == "list"){
				lapply(clust.branch, subcluster)
				}
	}	
			
		assign.mod.names <- function(clust.branch){
			if(class(clust.branch) != "list"){
				clust.name <- paste(clust.branch, collapse = sep)
				return(clust.name)
				}
			if(class(clust.branch) == "list"){
				lapply(clust.branch, assign.mod.names)
				}
		}	
			
	#==========================================
	
		
	final.clust <- make.clust.tree(clust)
	all.mod.size <- unlist(get.mod.size(final.clust))

	while(any(all.mod.size > max.mod.size)){
		final.clust <- subcluster(final.clust)
		all.mod.size <- unlist(get.mod.size(final.clust))
		}
	
	
	mod.name.tree <- assign.mod.names(final.clust)
	mod.names <- unlist(mod.name.tree)
	mod.num <- 1:length(mod.names)

	mod.pos <- sapply(mod.names, function(x) grep(x, mod.name.tree))
	clust.name.tree <- unlist(mod.name.tree, recursive = FALSE)
	while(class(clust.name.tree) == "list"){
		mod.pos <- cbind(mod.pos, sapply(mod.names, function(x) grep(x, clust.name.tree)))
		clust.name.tree <- unlist(clust.name.tree, recursive = FALSE)	
		}

	#pheatmap(mod.pos, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)

	if(class(mod.pos) == "integer"){
		mod.pos.mat <- matrix(mod.pos, ncol = 1)
		rownames(mod.pos.mat) <- names(mod.pos)
		mod.pos <- mod.pos.mat
	}
	colnames(mod.pos) <- paste("tier", 1:ncol(mod.pos), sep = "")


	#identify small module and group them together
	#module appropriately
	num.per.mod <- sapply(rownames(mod.pos), function(x) length(strsplit(x, "-")[[1]]))
	small.mods <- which(num.per.mod < min.mod.size)
	big.mods <- which(num.per.mod >= min.mod.size)
	if(length(small.mods) > 0){
		small.mod.mat <- mod.pos[small.mods,,drop=FALSE]
		big.mod.mat <- mod.pos[big.mods,,drop=FALSE]
		#pheatmap(small.mod.mat, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)
		#pheatmap(big.mod.mat, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)
		#collapse the small modules into a single module
		#it's okay if it's bigger than the maximum.
		all.small.genes <- unique(as.vector(unlist(sapply(rownames(small.mod.mat), function(x) strsplit(x, "-")[[1]]))))
		extra.mod.name <- paste(all.small.genes, collapse = sep)
		extra.row <- colMeans(small.mod.mat)
		final.mods <- rbind(big.mod.mat, extra.row)
		rownames(final.mods)[nrow(final.mods)] <- extra.mod.name
		#pheatmap(final.mods, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)

	}else{
		final.mods <- mod.pos
	}


	return(final.mods)

	}





