#This function is based on the code of mediation.analysis.
#It performs a mediation analysis.
#It asks for the variance of the outcome variable
#(out.var) explained by the explanatory variable 
#(exp.var) using a linear model (lm). It then places 
#each of the mediators in the linear model, and 
#gets the new variance explained by the exp.var
#after mediating on each mediator.
#out.var is a matrix of n samples in rows and
#v columns corresponding to different variables,
#for example fear conditioning, amyloid, etc.
#exp.var is a vector of length n samples. 
#mediators is a matrix with n rows and 
#m mediators
#Unlike mediation.analysis(), this function 
#takes covariates into account. It also performs
#the reverse causality tests. So, in addition to 
#putting the mediator as an explanatory variable in
#the model for the outcome variable, it also looks
#for explanation of the mediator by the outcome variable.
#The results of this function are reported differently
#than mediation.analysis(). In mediation.analysis()
#the results are reported as the remaining percent variance
#explained after the mediator is taken into effects.
#This function reports the difference between the 
#percent variance explained before and after the 
#mediator is accounted for.


#library(rgl)
#exp.var <- rnorm(100)
#out.var <- matrix(rnorm(300, sd = 0.7), nrow = 100, ncol = 3)
#out.var <- apply(out.var, 2, function(x) (exp.var*runif(1,0,2))+x)
#colnames(out.var) <- paste0("trait", 1:ncol(out.var))
#mediators <- matrix(rnorm(10000), nrow = 100, ncol = 100)
# mediators <- apply(mediators, 2, function(x) (exp.var*out.var[,runif(1, 1, 3)]-x*0.5))
# colnames(mediators) <- paste0("mediator", 1:ncol(mediators))
#plot3d(exp.var, out.var[,1], mediators[,1])
#if calc.mediator.exp.var is TRUE, this function will also calculate the
#percent variance explained by each mediator alone without the explanatory variable.

mediation.analysis2 <- function(out.var, exp.var, mediators, covar = NULL,
plot.results = FALSE){
	
	#========================================================
	#adjust variables for covariates.
	#========================================================
    if(!is.null(covar)){
        cat("Adjusting variables for covariates...\n")
        adj.out <- adjust(out.var, covar)
        adj.med <- adjust(mediators, covar)
    }else{
        adj.out <- out.var
        adj.med <- mediators
    }

	#========================================================
	#do some quick checks for variance in the input matrices
	#========================================================
	
	test.var <- round(apply(adj.out, 2, function(x) var(x, na.rm = TRUE)), 2)
	if(any(test.var == 0)){
		stop("some outcome variables have 0 variance.")
		}
	
	test.var <- round(apply(adj.med, 2, function(x) var(x, na.rm = TRUE)), 2)
	if(any(test.var == 0)){
		zero.locale <- which(test.var == 0)
		adj.med <- adj.med[,-zero.locale]
		warning(paste("removing", length(zero.locale), "mediators with zero variance"))
		}
	
	test.var <- round(apply(adj.med, 1, function(x) var(x, na.rm = TRUE)), 4)
	na.locale <- which(is.na(test.var))
	zero.locale <- which(test.var == 0)
	if(length(na.locale) > 0 || length(zero.locale) > 0){
		remove.locale <- unique(c(na.locale, zero.locale))
		adj.med <- adj.med[-remove.locale,]
		adj.out <- adj.out[-remove.locale,]
		exp.var <- exp.var[-remove.locale]
		warning(paste("removing", length(remove.locale), "individuals with zero expression variance"))
		}

	
	#========================================================
	#find the variance of the outcome measure explained by each
	#explanatory variable
	#========================================================
		
	# boxplot(adj.out[,2]~exp.var)
	factor.test <- apply(adj.out, 2, function(x) lm(x~exp.var))
	orig.pct.exp <- unlist(lapply(factor.test, function(x) var.exp(x, 1)))
	
    #also find the variance of each transcript explained by genotype
    cat("Calculating variance of each mediator explained by explanatory variable...\n")
    rev.orig.model <- apply(adj.med, 2, function(x) lm(x~exp.var))
    rev.orig.var <- unlist(lapply(rev.orig.model, function(x) var.exp(x, 1)))

	#========================================================
	#put each mediator in the model and calculate the new
	#variance explained
	#========================================================	

    #initialize all the matrices we will use
	
    #difference in trait variance explained by genotype 
    #before and after accounting for expression.
    diff.trait.var <- matrix(NA, nrow = ncol(adj.med), ncol = ncol(adj.out))
	
    #difference in transcript variance explained by genotype
    #after accounting for the trait
    diff.med.var <- matrix(NA, nrow = ncol(adj.med), ncol = ncol(adj.out))

    colnames(diff.trait.var) <- colnames(diff.med.var) <- colnames(adj.out)
	rownames(diff.trait.var) <- rownames(diff.med.var) <- colnames(adj.med)

	
	for(i in 1:ncol(adj.out)){
		cat("calculating mediators of", colnames(adj.out)[i], "...\n")
		mediator.models <- apply(adj.med, 2, function(x) lm(adj.out[,i]~x+exp.var))
        final.var.exp <- unlist(lapply(mediator.models, function(x) var.exp(x, 2)))
		diff.trait.var[,i] <- orig.pct.exp[i] - final.var.exp #record the drop in variance explained
	
        #calculate the variance of the mediator explained by the trait.
		cat("calculating reverse causality...\n")
        #next calculate the transcript variance explained when the trait is accounted for
		rev.mediator.models <- apply(adj.med, 2, function(x) lm(x~adj.out[,i]+exp.var))
        rev.final.var.exp <- unlist(lapply(rev.mediator.models, function(x) var.exp(x, 2)))
		diff.med.var[,i] <- rev.orig.var - rev.final.var.exp #record the drop in variance explained

        if(plot.results){
            par(mfrow = c(1,2))
            plot(rev.orig.var, rev.final.var.exp, ylab = "Variance Explained After Accounting for Trait",
            xlab = "Variance Explained by Genotype Alone")
            abline(0,1)

            plot(diff.trait.var[,i], diff.med.var[,i], 
            xlab = "Degree of Mediation of Transcript on Trait", 
            ylab = "Degree of Mediation of Trait on Transcript")
            abline(0,1);abline(h = 0, v = 0)
            mtext(colnames(out.var)[i], outer = TRUE, side = 3, line = -2)
        }
    }
	
    all.var.exp <- rbind(orig.pct.exp, diff.trait.var)
	rev.var.exp <- cbind(rev.orig.var, diff.med.var)
    colnames(rev.var.exp)[1] <- "Variance_of_Mediator_Explained"
	
    final.result <- list("forward" = all.var.exp, "reverse" = rev.var.exp)
    return(final.result)
	
}