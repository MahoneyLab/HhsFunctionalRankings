#This function makes a plot with two different y axes
#The code is adapted from https://www.r-bloggers.com/r-single-plot-with-two-different-y-axes/
#ylab1 = "vector1"; ylab2 = "vector2"; xlab = ""; col = c("#7fc97f", "#beaed4")
#legend.pos = "topleft"; legend.lab = c("vector1", "vector2"); plot.type = c("l", "l")
#lty = c(1,1); lwd = c(1,1); pch = c(NA, NA); cex = c(1,1)

plot.two.y.axes <- function(x1, y1, x2, y2, ylab1 = "vector1", ylab2 = "vector2", 
xlab = "", col = c("#7fc97f", "#beaed4"), legend.pos = "topleft", 
plot.type = c("l", "l"), lty = c(1,1), lwd = c(1,1), pch = c(NA, NA), cex = c(1,1)){
	
	d <- data.frame(cbind(x1, y1, x2, y2))

	par(mar = c(5,5,2,5))
	with(d, plot(x1, y1, type = plot.type[1], col = col[1], ylab = ylab1, lty = lty[1], lwd = lwd[1], 
	cex = cex[1], pch = pch[1], xlab = xlab))
	
	par(new = T)
	with(d, plot(x2, y2, type = plot.type[2], axes = F, xlab = NA, ylab = NA, lty = lty[1], 
	lwd = lwd[1], cex = cex[1], pch = pch[1], col = col[2]))
	axis(side = 4)
	mtext(side = 4, line = 3, ylab2)
	
	legend(legend.pos, legend = c(ylab1, ylab2), lty = lty, pch = pch, col = c(col[1], col[2]), lwd = lwd)
	
	
}