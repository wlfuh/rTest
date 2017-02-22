# Globals
random_seed <- sample(1:1000, 1)
masterlist <- list()
info <- 
  c( "Each index of masterlist represents the corresponding state",
     "predcs_prob = probability matrix",
     "assignments = Hungarian assignments based on probability",
     "predcs = Predicted Chemical shift data",
     "merged_cs = Merged predcs with corresponding cs from reference",
     "assign_list = List of each assignment per iteration",
     "predcs_list = List of predcited chemical shift data modified by noise due to error")

# Functions
assign <- function(x, y){
  suppressPackageStartupMessages(require(clue))
  # assigns elements in y to elements in x 
  # Input -- x (vector): e.g., chemical shifts 
  # Input -- y (vector): e.g., chemical shifts which are to be mapped to those in x (note that x can be greater than y)
  if(length(x) < length(y)){
    temp <- x
    x <- y
    y <- temp
  }
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
  costmat <- abs(xmat-ymat)
  a <- solve_LSAP(costmat)
  
  a_mat <- matrix(0,nrow=length(y),ncol=length(x),byrow=T)
  for (i in 1:length(y)){
    j <- as.vector(a)[i]
    a_mat[i,j] <-1
  }
  return(list(a=a,a_mat=a_mat,costmat=costmat))
}

get_assignment_cost <- function(a, costmat){
  require(clue)
  # get cost of assignment
  # Input -- a (vector): linear assignment object returned by solve_LSAP
  # Input -- costmat (matrix): cost matrix
  cost = as.integer(0)
  for(i in 1:length(a))
    cost = cost + costmat[i, a[i]]
  return(cost)
}

add_noise <- function(x, scale=1){
  # add random normal noise based on error
  x$cs+rnorm(1, 0, scale*x$error)
}

get_assigned_cs <- function(predcs, iter=5000, scale=0.75, dumpname="", numstates=11, saveRes=FALSE, freq_output=5000, refdata=refdata, output=output){
  #print("I am here")
  # calculate assigned chemical shift values to data
  if(dumpname == "")
    dumpname <- paste(random_seed,"_",iter,"_iterations_",scale,"_scale",sep="")
  
  predcs <- predcs[order(predcs$resid, predcs$nucleus),]
  
  if(nrow(predcs)>nrow(refdata)){refdata <- rbind(refdata,rbind(refdata, refdata))}
  
  prob <- matrix(0, nrow=nrow(predcs), ncol=nrow(refdata))
  testdata <- refdata[,"cs"]
  predcs_list <- NULL
  assign_list <- NULL
  state <- unique(predcs$state)
  
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  for(i in 1:iter){
    setTxtProgressBar(pb, i)    
    predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=scale)
    predcs_mod <- predcs_mod[order(predcs_mod$resid, predcs_mod$nucleus),]
    
    a <- assign(testdata, predcs_mod$V1)
    prob = prob + a$a_mat
    if (i%%freq_output==0){
			predcs_prob <- prob / i
			a <- solve_LSAP(1 - predcs_prob)
	
			# sets the assignment chemical shift to predicted chemical shift
			predcs$assigned <- refdata[a,"cs"]	
			tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
			tmp <- tmp[order(tmp$resid, tmp$nucleus),]
			tmp <- rename(tmp, c("resname.x" = "resname", "cs.x"="actual", "cs.y"="predicted"))
			tmp <- unique(tmp[, c("state", "resid", "resname", "nucleus", "actual", "assigned", "predicted")])
			write.table(tmp, file = paste(output, "_", i, "_", state, ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
			write.table(predcs_prob, file = paste(output, "_prob_", i, "_", state, ".RData", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
		}    
  }
  return(TRUE)
  
	#predcs_prob <- prob / i
	#a <- solve_LSAP(1 - predcs_prob)
	## sets the assignment chemical shift to predicted chemical shift
	#predcs$assigned <- refdata[a,"cs"]	
	#tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
	#tmp <- tmp[order(tmp$resid, tmp$nucleus),]
	#tmp <- rename(tmp, c("resname.x" = "resname", "cs.x"="actual", "cs.y"="predicted"))
	#tmp <- unique(tmp[, c("state", "resid", "resname", "nucleus", "actual", "assigned", "predicted")])	
	#tmpout <- paste(output, "_", i, "_", state, ".txt", sep = "")
	#write.table(tmp, file = tmpout, quote = FALSE, row.names = FALSE, col.names = FALSE)
	#tmpout <- paste(output, "_prob_", i, "_", state, ".RData", sep = "")
	#save(file = tmpout, predcs_prob)  
  #close(pb)  
  #tmp
}

matchModel <- function(iter=5000, scale=0.75, method="mean",filename="data/predicted_CS_table_test_clean.txt",label="Test", saveRes=FALSE, freq_output = 5000, refdata=refdata, output=output, parallel = FALSE){
  modeldata <- read.table(filename, header=FALSE)  
  if (ncol(modeldata==5)){names <- c("state", "resid", "resname", "nucleus", "cs")}
  if (ncol(modeldata==6)){names <- c("state", "resid", "resname", "nucleus", "cs", "na")}       
  colnames(modeldata) <- names
  
  if(parallel){
  	modeldata <- modeldata[ 1:nprocessors %in% modeldata$state, ]
  }
  
  modeldata <- modeldata[,names]
  modeldata <- merge(modeldata, accu)
  levels <- length(unique(modeldata$state))
  modeldata['assigned'] <- NA
  ddply(.data=modeldata, .var=c("state"),.fun = get_assigned_cs, iter=iter, scale=scale, dumpname=label, saveRes=saveRes, freq_output=freq_output, , refdata=refdata, output=output, .parallel = parallel)  
}
