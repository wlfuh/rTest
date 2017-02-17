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
    temp = x
    x = y
    y = temp
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

get_assigned_cs <- function(predcs, iter=5000, scale=0.75, dumpname="", numstates=11, saveRes=FALSE){
  # calculate assigned chemical shift values to data
  if(dumpname == "")
    dumpname <- paste(random_seed,"_",iter,"_iterations_",scale,"_scale",sep="")
  
  predcs <- predcs[order(predcs$resid, predcs$nucleus),]
  prob <- matrix(0, nrow=nrow(predcs), ncol=nrow(refdata))
  testdata <- refdata[,"cs"]
  predcs_list <- NULL
  assign_list <- NULL
  
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  for(i in 1:iter){
    setTxtProgressBar(pb, i)    
    predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=scale)
    predcs_mod <- predcs_mod[order(predcs_mod$resid, predcs_mod$nucleus),]
    
    a <- assign(testdata, predcs_mod$V1)
    prob = prob + a$a_mat
    if(saveRes){
      predcs_list <- c(predcs_list, predcs_mod)
      assign_list <- c(assign_list, a)
    }
  }
  close(pb)
  
  predcs_prob <- prob / iter
  a <- solve_LSAP(1 - predcs_prob)
  
  # sets the assignment chemical shift to predicted chemical shift
  predcs$assigned <- refdata[a,"cs"]
  
  tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
  tmp <- tmp[order(tmp$resid, tmp$nucleus),]
  #print(tmp)
  if(saveRes){
    masterlist[[predcs$state[1]]] <<- list(predcs_prob=predcs_prob, assignments=a, predcs=predcs, merged_cs=tmp,
                                           assign_list=assign_list, predcs_list=predcs_list)
    if(predcs$state[1] == numstates){
      save_data(paste(dumpname,"_",iter,"_iterations_",scale,"_scale_data",sep=""), folder=dumpname)
    }
  }
  tmp <- rename(tmp, c("resname.x" = "resname", "cs.x"="actual", "cs.y"="cs"))
  tmp[, c("state", "resid", "resname", "nucleus", "cs", "assigned", "actual")]
}

save_data <- function(filename, folder=""){
  if(folder != ""){
    dir.create(paste("output/",folder,sep=""), showWarnings = FALSE)
    filename <- paste(folder,"/" ,filename, sep="")
  }
  save(info, masterlist, file=paste("output/",filename,".RData",sep=""))
}

matchModel <- function(iter=5000, scale=0.75,method="mean",filename="data/predicted_CS_table_test_clean.txt",label="Test", saveRes=FALSE){
  modeldata <- read.table(filename, header=FALSE)
  if (ncol(modeldata==5)){names <- c("state", "resid", "resname", "nucleus", "cs")}
  if (ncol(modeldata==6)){names <- c("state", "resid", "resname", "nucleus", "cs", "na")}
  colnames(modeldata) <- names
  modeldata <- modeldata[,names]
  modeldata <- merge(modeldata, accu)
  levels <- length(unique(modeldata$state))
  modeldata['assigned'] <- NA
  temp <- ddply(.data=modeldata, .var=c("state"),.fun = get_assigned_cs, iter=iter, scale=scale, dumpname=label, saveRes=saveRes)
  return(temp)
}
