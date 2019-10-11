#This function takes x alone as a list, or x with a factor y 
#generates a boxplot with a p value from an ANOVA.

plot.box.with.model <- function(x, y = NULL, xlim = NULL, ylim = NULL, col = "black", 
pch = 16, main = "", group.names = NULL, ylab = "Y", add.stripchart = TRUE){
	
	if(class(x) == "list"){
		factors <- 1:length(x)
		y.list <- mapply(function(a,b) rep(a, length(b)), factors, x)
		y <- unlist(y.list)
		x <- unlist(x)
	}
	

	if(is.null(group.names)){group.names = levels(as.factor(y))}
	model <- aov(x~y)
	model.vals <- anova(model)
	p <- signif(model.vals[1,ncol(model.vals)], 3)
	boxplot(x~as.factor(y), main = paste(main, "\np =", p), ylab = ylab, names = group.names, 
	xlim = xlim, ylim = ylim)
	
	if(add.stripchart){
		y.vals <- unique(y)
		val.pos <- 1:length(y.vals)
		y.val.pos <- rep(NA, length(y))
		for(i in 1:length(val.pos)){
			y.ind <- which(y == y.vals[i])
			y.val.pos[y.ind] <- rep(val.pos[i], length(y.ind))
		}
		points(jitter(y.val.pos),x, pch = pch, col = col)
	}	

invisible(p)		
}