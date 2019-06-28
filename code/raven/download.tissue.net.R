#I'm still working out how best to do this function
#we want to be able to 

download.tissue.net <- function(tissue = NULL, organism = c("mouse", "human"), top.edges.only = TRUE, project.dir = "."){

	require(XML)
	require(RCurl)
		
	cur.dir <- getwd()

	if(length(organism) == 2){
		organism <- "mouse"	
		}
		
	if(is.null(tissue)){
	if(organism == "mouse"){
		txt <- getURL("http://fntm.princeton.edu/download/")
		}else{
		txt <- getURL("http://hb.flatironinstitute.org/download/")	
		}
		
	tree <- htmlTreeParse(txt, asText = TRUE)
	tissue.table <- tree$children$html[[2]][[2]][[3]][[2]]

	tissues <- rep(NA, length(tissue.table))
	for(i in 1:length(tissue.table)){
		tissue.v <- as.character(tissue.table[[i]][[1]][[1]])
		tissues[i] <- tail(tissue.v, 1)
		}
	
	cat("Please select from the following tissues...\n")
	cat(tissues, sep = "\n")
	}else{
	
	if(top.edges.only){
		tissue.file <- file.path(project.dir, paste0(gsub(" ", "_", tissue), "_top.gz", sep = ""))
		}else{
		tissue.file <- paste0(gsub(" ", "_", tissue), ".gz", sep = "")	
		}

	if(organism == "mouse"){
		base.url <- "http://fntm.princeton.edu/static//networks10"
		}else{
		base.url <- "https://s3-us-west-2.amazonaws.com/humanbase/networks/"
		}
		
	download.file(paste(base.url, tissue.file, sep = "/"), destfile = tissue.file)
	
	cat("Unzipping network file...\n")
	system(paste("gunzip", tissue.file))
	cat("Reading in network file...\n")
	tissue.net <- read.table(gsub(".gz", "", tissue.file))
	cat("Saving R binary version of network file...\n")
	saveRDS(tissue.net, gsub(".gz", ".RData", tissue.file))
	unlink(gsub(".gz", "", tissue.file))

	invisible(tissue.net)
	setwd(cur.dir)
	}
}