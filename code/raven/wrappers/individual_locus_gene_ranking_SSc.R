#this is a wrapper for ranking genes in individual loci 
library(biomaRt)
library(stringr)
library(vioplot)
mart = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
snp.db = useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

cape.dir <- "~/Documents/git_repositories/capeDO"
setwd(cape.dir)
all.fun <- list.files(pattern = "*.R")
for(i in 1:length(all.fun)){source(all.fun[i])}

useful.dir <- "~/Documents/git_repositories/useful_r_code"
setwd(useful.dir)
all.fun <- list.files(pattern = "*.R")
for(i in 1:length(all.fun)){source(all.fun[i])}

code.dir <- "~/Documents/git_repositories/raven"
setwd(code.dir)
all.fun <- list.files(pattern = "*.R")
for(i in 1:length(all.fun)){source(all.fun[i])}

data.dir <- "~/Documents/Grants/R21/Scleroderma/Results/capeRel_testing"
setwd(data.dir)
cape.results <- readRDS("pop.RData")
geno.obj <- readRDS("popGeno.RData")

project.dir <- "~/Documents/Grants/R21/Scleroderma/Results/RAVEN"
results.dir <- "~/Documents/Grants/R21/Scleroderma/Results/RAVEN/networks"
tissue.name = "skin"
organism = "human"
phenotype.name = "antibody"
upstream.buffer = 2.5e6
downstream.buffer = 2.5e6

setwd(project.dir)
# tissue.net  <- download.tissue.net(tissue.name, organism, TRUE)
tissue.net <- readRDS(paste(tissue.name, "_top.RData", sep = ""))

int.loci <- as.matrix(read.table("Motifs.of.Interest.txt", sep = "\t", header = TRUE))

mp.genes <- get.mp.genes(mp.term = phenotype.name, organism)

#====================================================================
# build the tissue-specific adjacency matrix for phenotype-associated 
# genes
#====================================================================

tissue.mat <- tissue.adj.mat(tissue.net = tissue.net, gene.list = mp.genes[[3]], inc.all.genes = TRUE)
# saveRDS(tissue.mat, paste(tissue.name, "Adjacency.RData", sep = ""))
# tissue.mat <- readRDS(paste(tissue.name, "Adjacency.RData", sep = ""))
tissue.mat <- remove.unconnected.tp(tissue.mat = tissue.mat)

plot.tp.adj(tissue.mat)
tp.mat <- get.tp(tissue.mat)
hist(tp.mat[which(tp.mat != 0)])
pdf(paste0("tp.mat.", tissue.name, ".pdf"), height = 10, width = 10)
heatmap(tp.mat)
dev.off()


tp.mat <- tissue.mat[,as.character(rownames(tissue.mat))]
neg.genes <- setdiff(colnames(tissue.mat), rownames(tissue.mat))
tn.mat <- tissue.mat[,as.character(neg.genes)]

num.tp <- ncol(tp.mat)
num.tn <- ncol(tn.mat)

full.data <- cbind(tp.mat, tn.mat)
labels <- as.factor(c(rep("pos", ncol(tp.mat)), rep("neg", ncol(tn.mat))))


orig.cost.list <- 10^seq(-5, 3, 1)
tuned.cost.list <- tune.cost.parameter(tissue.mat, n.trials = 5, C.list = orig.cost.list, verbose = TRUE)

ntrials <- 100
svm.pheno.genes <- tissue.svm(tissue.mat = tissue.mat, n.trials = ntrials, C.list = tuned.cost.list, vote.min = 0.5, verbose = TRUE)
saveRDS(svm.pheno.genes, paste0("svmPhenoGenes.", tissue.name, ".RData"))
# svm.pheno.genes <- readRDS(paste0("svmPhenoGenes.", tissue.name, ".RData"))

if(!file.exists("Sampled.Predictions.RData")){
	all.predict <- matrix(NA, nrow = ntrials, ncol = ncol(full.data))
	sample.labels <- as.factor(c(rep("pos", num.tp), rep("neg", num.tp)))
	
	for(n in 1:ntrials){
		print(n)
		sample.mat <- t(cbind(tp.mat, tn.mat[,sample(1:num.tn, num.tp)]))
		results <- cv.linear.svm(data.mat = sample.mat, data.labels = sample.labels, C.list = 4^seq(-5, -3, 1))
		prediction <- predict(results$opt.model, newdata = t(full.data))
		all.predict[n,] <- prediction
		}
	saveRDS(all.predict, "Sampled.Predictions.RData")
	}else{
	predictions <- readRDS("Sampled.Predictions.RData")
	}
gene.votes <- apply(predictions, 2, function(x) length(which(x == "pos")))

vote.order <- order(gene.votes, decreasing = TRUE)
ordered.votes <- gene.votes[vote.order]
plot(ordered.votes, type = "p", pch = 16, cex = 0.2)

#====================================================================
# where do the locus genes fall in the voting list?
#====================================================================
get.id.vote <- function(entrez.id, ordered.list){
	if(is.na(entrez.id)){return(NA)}
	locale <- which(names(ordered.list) == entrez.id)
	if(length(locale) > 0){return(locale)}else{return(length(ordered.list))}
	}

gene.expr.cor <- function(gene1, gene2){
	gene1.idx <- which(gene.names[,2] == str_to_upper(gene1))
	gene2.idx <- which(gene.names[,2] == str_to_upper(gene2))
	
	if(length(gene1.idx) > 0 && length(gene2.idx) > 0){

	gene1.id <- gene.names[gene1.idx,1]
	gene2.id <- gene.names[gene2.idx,1]
	
	gene1.id.locale <- which(entrez.id == gene1.id)
	gene2.id.locale <- which(entrez.id == gene2.id)
	
	gene1.expr <- apply(just.expr[,gene1.id.locale,drop=FALSE], 1, mean)
	gene2.expr <- apply(just.expr[,gene2.id.locale,drop=FALSE], 1, mean)
	

	model <- lm(gene2.expr~gene1.expr)
	model.sum <- coef(summary(model))
	plot(gene1.expr, gene2.expr, main = paste(gene1, gene2, "\nr =", signif(cor(gene1.expr,gene2.expr), 2), "\np =", signif(model.sum[2,4], 2)), col = ssc.col, pch = 16, xlab = gene1, ylab = gene2)
	abline(model)	
	par(xpd = TRUE)
	plot.height <- max(gene2.expr, na.rm = TRUE) - min(gene2.expr, na.rm = TRUE) 
	plot.width <- max(gene1.expr, na.rm = TRUE) - min(gene1.expr, na.rm = TRUE) 
	legend.x <- min(gene1.expr, na.rm = TRUE) - plot.width*0.13
	legend.y <- min(gene2.expr, na.rm = TRUE) - (plot.height*0.18)
	legend(legend.x, legend.y, legend = ssc.col.names, col = ssc.cols, pch = 16, horiz = TRUE)
	par(xpd = FALSE)
	return(c("r" = cor(gene1.expr,gene2.expr), "p" = model.sum[2,4]))
	}else{
		return(c("r" = NA, "p" = NA))
		}
	
	}

get.gene.pairs <- function(good.pair.list.element){
	
	genes1 <- good.pair.list.element[[1]][,1]
	genes2 <- good.pair.list.element[[2]][,1]
	num.genes2 <- length(genes2)
	pair.mat <- matrix(NA, nrow = length(genes1)*length(genes2), ncol = 2)
	pair.mat[,2]  <- rep(genes2, length(genes1))
	pair.mat[,1] <- rep(genes1, each = length(genes2))
	return(pair.mat)
	}

plot.ssc.gene.expr <- function(gene.id = NULL, gene.name = NULL, plot.type = c("box", "violin")){
	plot.type <- plot.type[1]
	if(is.null(gene.id)){
		gene.id <- getBM("entrezgene", "external_gene_name", gene.name, mart)[1,1]
		}
	if(is.null(gene.name)){
		gene.name <- getBM("external_gene_name", "entrezgene", gene.id, mart)[1,1]
		}
	gene.locale <- which(entrez.id == gene.id)
	if(length(gene.locale) > 0){
		grp.test <- anova(aov(just.expr[,gene.locale]~grp))
		ssc.test <- t.test(just.expr[,gene.locale]~sc.v.n)
		pvals <- grp.test$"Pr(>F)"
		pval <- tail(pvals[which(!is.na(pvals))], 1)
		all.p[i] <- ssc.test$p.value
		par(mfrow = c(1,2))

		u_grp <- unique(grp)
		grp.list <- lapply(u_grp, function(x) just.expr[which(grp == x),gene.locale])
		names(grp.list) <- u_grp
		if(plot.type == "box"){
			boxplot(grp.list, main = paste(gene.name, "\np =", signif(pval, 2)), las = 2)
			}else{
			vioplot(grp.list[[1]], grp.list[[2]], grp.list[[3]], grp.list[[4]], grp.list[[5]], col = "gray", names = u_grp)
			mtext(paste(gene.name, "\np =", signif(pval, 2)))
			stripchart(grp.list, vertical = TRUE, pch = 16, add = TRUE, method = "jitter", col = "palevioletred1")
			}
		
		u_sc.n <- unique(sc.v.n)
		sc.n.list <- lapply(u_sc.n, function(x) just.expr[which(sc.v.n == x),gene.locale])
		names(sc.n.list) <- u_sc.n
		if(plot.type == "box"){
			boxplot(sc.n.list, main = paste(gene.name, "\np =", signif(ssc.test$p.value, 2)))
			}else{
			vioplot(sc.n.list[[1]], sc.n.list[[2]], col = "gray", names = u_sc.n)
			mtext(paste(gene.name, "\np =", signif(ssc.test$p.value, 2)))
			stripchart(sc.n.list, vertical = TRUE, pch = 16, add = TRUE, method = "jitter", col = "palevioletred1")
			}

		}else{
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		text(0.5, 0.5, paste("No expression data for", gene.name))
		}
	}

gene.expr <- function(gene.id){
	gene.locale <- which(entrez.id == gene.id)
	ssc.locale <- which(sc.v.n == "SSc")
	cont.locale <- which(sc.v.n == "Normal")
	if(length(gene.locale) > 0){
		ssc.vals <-just.expr[ssc.locale,gene.locale]
		cont.vals <-just.expr[cont.locale,gene.locale]
		return(c("SSc" = mean(ssc.vals), "Cont" = mean(cont.vals)))		
		}else{
		return(c("SSc" = NA, "Cont" = NA))		
		}	
	}

# hand.picked <- read.table("~/Desktop/candidate_pairs_by_cape_plot.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hand.picked <- int.loci
u_snps <- unique(c(hand.picked[,1], hand.picked[,2]))
for(i in 1:length(u_snps)){
	print(i)
	snp <- u_snps[i]
	nearby.genes <- gene.near.snp(snp.name = snp, mart, snp.db, upstream.buffer, downstream.buffer)
	not.na.locale <- which(!is.na(nearby.genes[[2]][,1]))
	nearby.genes <- nearby.genes[[2]][not.na.locale,]
	gene.ids <- nearby.genes[,1]
	#nearby.genes[,2]
	if(length(gene.ids) > 0 && gene.ids[1] != "cannot find SNP"){
		gene.rank <- unlist(lapply(gene.ids, function(x) get.id.vote(x, ordered.votes)))
		gene.vote <- ordered.votes[gene.rank]
		rank.table <- cbind(gene.rank, gene.vote, nearby.genes)
		rank.table <- rank.table[order(rank.table[,1]),]
		ranks <- condense.table(tableX = rank.table, condense.by = 3, col.to.collapse = c(1:2,4:7), col.to.concat = 9)
		write.table(ranks, paste0("Ranked.Genes.", snp, ".", tissue.name, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)	
		}
	}


#expanded network
candidate.genes <- get.candidate.genes(90, 1)
u_candidate.genes <- unique(Reduce(rbind, lapply(candidate.genes, function(x) x[,1:3])))
expanded.tissue.mat <- tissue.adj.mat(tissue.net = tissue.net, gene.list = c(u_candidate.genes[,2], mp.genes[[3]]), inc.all.genes = TRUE)


#look at tissue-specific subnetworks of each gene
#how many of their neihbors are phenotype specific
overlap.with.pheno.genes <- vector(mode = "list", length = nrow(u_candidate.genes))
names(overlap.with.pheno.genes) <- u_candidate.genes[,1]

for(i in 1:nrow(u_candidate.genes)){
	print(i)
	subnet <- gene.subnetwork(entrezid = as.numeric(u_candidate.genes[i,2]), tissue.adj.mat = expanded.tissue.mat, num.steps = 2)

	# pdf(paste0("Gene.Subnetwork.", u_candidate.genes[i,1], ".pdf"), width = 10, height = 10)
	# v.col <- rep("lightblue", vcount(subnet))
	# v.col[which(V(subnet)$name %in% mp.genes[[2]])] <- "orange"
	# v.col[which(V(subnet)$name == u_candidate.genes[i,1])] <- "red"
	# plot(subnet, vertex.size = 5, layout = layout.kamada.kawai, edge.width = E(subnet)$weight*2, vertex.color = v.col, vertex.label.dist = 0.3)
	# hist(E(subnet)$weight, breaks = 50, main = "Histogram of Edge Weights")
	
	# enrichment <- gprofiler(V(subnet)$name, "hsapiens")
	# plot.enrichment(enrichment, text.size = 1)
	# dev.off()
	if(is.igraph(subnet)){
		overlap.with.pheno.genes[[i]] <- intersect(V(subnet)$name, mp.genes[[2]])
		}
	
	}
	
num.overlap <- unlist(lapply(overlap.with.pheno.genes, length))
	
#look at expression differences of candidate genes
expr <- read.table("~/Documents/Data/Scleroderma/data_files/Expression/Hinchcliff.lpcl", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
expr <- expr[-1,]
oligo.data <- expr[,1:2]
split.id <- strsplit(oligo.data[,2], "\\^")
entrez.id <- trimws(unlist(lapply(split.id, function(x) x[length(x)])))
gene.names <- getBM(c("entrezgene", "external_gene_name"), "entrezgene", entrez.id, mart)

just.expr <- t(expr[,4:ncol(expr)])
grps <- read.table("~/Documents/Data/Scleroderma/data_files/Expression/Hinchcliff.phen", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
grp <- grps[2:nrow(grps),2]
sc.v.n <- grp
sc.v.n[which(sc.v.n != "Normal")] <- "SSc"

all.p <- rep(NA, nrow(u_candidate.genes))
pdf("Candidate.Gene.Expression.Handpicked.pdf", width = 10, height = 5)
for(i in 1:nrow(u_candidate.genes)){
	plot.ssc.gene.expr(gene.id = u_candidate.genes[i,2], gene.name = u_candidate.genes[i,1], plot.type = "box")
	}
dev.off()

#color by SSc and Control
ssc.cols <- c("black", "blue")
ssc.col.names <- unique(sc.v.n)
ssc.col.labels <- sc.v.n


#color by group
ssc.cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")
ssc.col.names <- unique(grp)
ssc.col.labels <- grp

ssc.col <- rep(NA, length(sc.v.n))
for(i in 1:length(ssc.col.names)){
	ssc.col[which(ssc.col.labels == ssc.col.names[i])] <- ssc.cols[i]
	}

candidate.table <- cbind(u_candidate.genes, all.p)
write.table(candidate.table, "Candidate.Info.Handpicked.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#find SNP pairs with highly ranked nearby genes
good.pairs <- 1
good.pair.list <- list()
vote.cutoff = 90
for(i in 1:nrow(int.loci)){
	snp1.locale <- which(names(candidate.genes) == int.loci[i,1])
	snp2.locale <- which(names(candidate.genes) == int.loci[i,2])
	
	if(length(snp1.locale) > 0 && length(snp2.locale) > 0){
		cand.gene1 <- candidate.genes[[snp1.locale]]
		cand.gene2 <- candidate.genes[[snp2.locale]]
		
		above.cutoff1 <- which(cand.gene1[,"gene.vote"] > vote.cutoff)
		above.cutoff2 <- which(cand.gene2[,"gene.vote"] > vote.cutoff)
		
		if(length(above.cutoff1) > 0 && length(above.cutoff2) > 0){
			cand1 <- cand.gene1[above.cutoff1,,drop=FALSE]
			cand1.genes.locale <- which(candidate.table[,1] %in% cand1[,1])
			cand1 <- cbind(cand1, candidate.table[cand1.genes.locale,4])
			
			cand2 <- cand.gene2[above.cutoff2,,drop=FALSE]
			cand2.genes.locale <- which(candidate.table[,1] %in% cand2[,1])
			cand2 <- cbind(cand2, candidate.table[cand2.genes.locale,4])


			cand.pair <- list(cand1, cand2)	
			names(cand.pair) <- int.loci[i,1:2]
			good.pair.list[[good.pairs]] <- cand.pair
			names(good.pair.list)[good.pairs] <- paste(int.loci[i,], collapse = " : ")
			good.pairs <- good.pairs+1
			}
		}
	}
	
	length(good.pair.list)
	saveRDS(good.pair.list, "good.pair.list.RData")
	
	# more.locale <- grep("more_than_additive", names(good.pair.list))
	
	# i <- 4
	# i = i + 1
	# good.pair.list[[more.locale[i]]]
source('~/Documents/R/shiny/cape_interaction_viewer/code/plot.bars.R')
pdf("All.motif.bars.pdf", width = 10, height = 5)
	par(mfrow = c(1,3), mar = c(5,5,10,5))
	for(i in 1:length(good.pair.list)){
		snp.names <- names(good.pair.list[[i]])
		snp.locale <- which(colnames(geno.obj$geno) %in% snp.names)
		marker.geno <- geno.obj$geno[,snp.locale]
		binned.geno <- marker.geno
		binned.geno[which(binned.geno > 0)] <- 1
		for(p in 1:ncol(cape.results$pheno)){
			plot.bars(binned.geno, snp.names, cape.results$pheno[,p], pheno.name = colnames(cape.results$pheno)[p], ref.centered = TRUE, error.bars = "se")
			}
		mtext(paste(snp.names, collapse = ", "), line = -2, outer = TRUE)
	}
	dev.off()
	
	possible.pairs <- unique(Reduce(rbind, lapply(good.pair.list, get.gene.pairs)))
	possible.pairs <- possible.pairs[which(apply(possible.pairs, 1, function(x) x[1] != x[2])),]
	
	all.cor <- matrix(NA, nrow = nrow(possible.pairs), ncol = 2)
	colnames(all.cor) <- c("r", "p")
	pdf("All.Candidate.Gene.Coexpression.Handpicked.pdf", width = 5, height = 5)
	for(i in 1:nrow(possible.pairs)){
		all.cor[i,] <- gene.expr.cor(gene1 = possible.pairs[i,1], gene2 = possible.pairs[i,2])
		}
	dev.off()
	
	gene.cor.table <- cbind(possible.pairs, all.cor)
	write.table(gene.cor.table, "Expression.Correlation.Table.Handpicked.txt", sep = "\t" ,quote = FALSE, row.names = FALSE)

#individual gene plotting by hand
gene.name <- c("SLC43A1", "GPSM3", "ARGLU1")
gene.ids <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", gene.name, mart)
pdf("Test.Gene.Expr.pdf", width = 10, height = 5)
for(i in 1:nrow(gene.ids)){
	if(!is.na(gene.ids[i,2])){
		plot.ssc.gene.expr(gene.id = gene.ids[i,2], gene.name = gene.ids[i,1])	
		}
	}
dev.off()


gene.expr.cor(gene1 = "ARGLU1", gene2 = "GPSM3")


#also look at genes nearest each SNP
near.source <- vector(mode = "list", length = nrow(int.loci))
near.target <- vector(mode = "list", length = nrow(int.loci))
names(near.source) <- names(near.target) <- int.loci[,1]
for(i in 1:nrow(int.loci)){
	report.progress(i, nrow(int.loci))
	near.source[[i]] <- gene.near.snp(snp.name = int.loci[i,1], mart, snp.db, 1e6, 1e6, TRUE, TRUE)
	near.target[[i]] <- gene.near.snp(snp.name = int.loci[i,2], mart, snp.db, 1e6, 1e6, TRUE, TRUE)
	}

saveRDS(near.source, "Genes.Near.Source.SNPs.RData")
saveRDS(near.target, "Genes.Near.Target.SNPs.RData")

source.genes <- unlist(lapply(near.source, function(x) if(length(x) == 2){x[[2]][2]}else{NA}))
target.genes <- unlist(lapply(near.target, function(x) if(length(x) == 2){x[[2]][2]}else{NA}))

nearest.edgelist <- cbind(source.genes, target.genes)
write.table(nearest.edgelist, "Edgelist.Nearest.Gene.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

nearest.edgelist <- read.table("Edgelist.Nearest.Gene.txt", stringsAsFactors = FALSE)
edge.cols <- rep("gray", nrow(nearest.edgelist))
edge.cols[which(int.loci[,"Class"] == "more_than_additive")] <- "orange"
edge.cols[which(int.loci[,"Class"] == "extreme")] <- "red"
edge.cols[which(int.loci[,"Class"] == "less_than_additive")] <- "blue"


na.locale <- which(is.na(nearest.edgelist), arr.ind = TRUE)
if(nrow(na.locale) > 0){
	nearest.edgelist <- nearest.edgelist[-na.locale,]
	}

edge.cols <- edge.cols[-na.locale]


net <- graph_from_edgelist(as.matrix(nearest.edgelist))
E(net)$color <- edge.cols
deg <- degree(net)

plot(net, vertex.size = sqrt(deg)*2, vertex.label.cex = 0.5, edge.arrow.size = 0.3)
