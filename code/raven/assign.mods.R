#This function takes the output of iter.cluster2 and the network that was
#clustered and retuns a vector assigning each node of the network to one
#of the clusters.

assign.mods <- function(net, mods, tier = ncol(mods)){

    nodes <- V(net)$name
    mod.v <- as.vector(mods[,tier])
    u_cluster <- sort(unique(mod.v))
    mod.assig <- rep(NA, length(nodes))

    for(i in 1:length(u_cluster)){
        cluster.locale <- which(mod.v == u_cluster[i])
        cluster.nodes <- rownames(mods)[cluster.locale]
        split.nodes <- as.vector(unlist(sapply(cluster.nodes, function(x) strsplit(x, "-"))))
        node.locale <- match(split.nodes, nodes)
        mod.assig[node.locale] <- u_cluster[i]
    }

    return(mod.assig)
}