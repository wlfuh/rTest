library(pracma)
library(plyr)
source("scripts/assignments.R")

# read in reference data
refdata <- read.table("data/measured_shifts_2KOC.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","dummy")

# add random normal noise based on error
add_noise <- function(x, scale=1){
  x$cs+rnorm(1, 0, scale*x$error)
}

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

# calculate assigned chemical shift values to data
get_assigned_cs <- function(predcs, iter=5000, scale=0.75, dumpname="", numstates=11){
  if(dumpname == "")
    dumpname <- paste(random_seed,"_",iter,"_iterations_",scale,"_scale",sep="")
  
  #dir.create(paste("output/",dumpname,sep=""), showWarnings = FALSE)
  #dumpname <- paste("output/",dumpname,"/state_",predcs$state[1],".RData",sep="")
  
  predcs <- predcs[order(predcs$resid, predcs$nucleus),]
  
  prob <- matrix(0, nrow=nrow(predcs), ncol=nrow(refdata))
  testdata <- refdata[,"cs"]
  
  predcs_list <- NULL
  assign_list <- NULL
  for(i in 1:iter){
    predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=scale)
    predcs_mod <- predcs_mod[order(predcs_mod$resid, predcs_mod$nucleus),]
    
    a <- assign(testdata, predcs_mod$V1)
    
    prob = prob + a$a_mat
    predcs_list <- c(predcs_list, predcs_mod)
    assign_list <- c(assign_list, a)
  }
  
  predcs_prob <- prob / iter
  a <- solve_LSAP(1 - predcs_prob)
  
  # sets the assignment chemical shift to predicted chemical shift
  predcs$assigned <- refdata[a,"cs"]
  
  tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
  tmp <- tmp[order(tmp$resid, tmp$nucleus),]
  
  # TODO, make all the states save in one file, idea: keep global variable that store each result in a list
  #save(predcs_prob, a, predcs, tmp, assign_list, predcs_list, file=dumpname)
  
  masterlist[[predcs$state[1]]] <<- list(predcs_prob=predcs_prob, assignments=a, predcs=predcs, merged_cs=tmp,
                                        assign_list=assign_list, predcs_list=predcs_list)
  if(predcs$state[1] == numstates){
    save_data(paste(dumpname,"_",iter,"_iterations_",scale,"_scale_data",sep=""), folder=dumpname)
  }
  out <- tmp[,c(7,1,3,2,4,10,8)]
  rename(out, c("resname.x" = "resname", "cs.x"="cs", "cs.y"="actual"))
}

save_data <- function(filename, folder=""){
  if(folder != ""){
    dir.create(paste("output/",folder,sep=""), showWarnings = FALSE)
    filename <- paste(folder,"/" ,filename, sep="")
  }
  save(info, masterlist, file=paste("output/",filename,".RData",sep=""))
}


# Perform iteration of assignments and returns most probable assignment
hung_trials <- function(predcs, iter=5000, scale=0.75,method="mean"){
  prob <- matrix(0, nrow=nrow(predcs), ncol=nrow(refdata))
  testdata <- refdata[,"cs"]
  
  for(i in 1:iter){
    predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=scale)
    predcs_mod <- predcs_mod[order(predcs_mod$resid, predcs_mod$nucleus),]
    
    a <- assign(testdata, predcs_mod$V1)
    
    prob = prob + a$a_mat
  }
  
  predcs_prob <- prob / iter
  a <- solve_LSAP(1 - predcs_prob)
  
  # sets the assignment chemical shift to predicted chemical shift
  predcs$assigned <- refdata[a,"cs"]
  
  # determine actual assignment for given residue id and nucleus
  tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
  
  func <- NULL
  if(method == 'mean'){
    func <- function(x){mean(abs(x$cs.x-x$assigned))}
  }
  else if(method == 'cor'){
    func <- function(x){cor(x$cs.x,x$assigned)}
  }
  else if(method == 'stderr'){
    func <- function(x){mean(abs((x$cs.x-x$assigned)/ x$error))}
  }
  else{
    print("Invalid method")
    return(NULL)
  }
  abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = func)
  
  return(get_correlation(list(cs=predcs, assigned=abc, iter=iter, scale=scale, assign_error=mean(abc$V1))))
}

# get correlation data between assigned and actual
get_correlation <- function(res, nuclei=c("H","C", "")){
  
  tmp <- merge(res$cs,refdata, by=c("resid","nucleus"),suffixes=c(".pred",".real"))
  dataCorr <- data.frame(nuc="B1",r=0,rsquare=0,RMSE=0,MAE=0,stringsAsFactors=FALSE)
  subtmp <- NULL
  for(nuc in nuclei){
    if(nuc == "")
      subtmp <- tmp
    else
      subtmp <- tmp[substr(tmp$nucleus,1,1)==nuc,]
    modeltmp <- lm(assigned~cs.real, data=subtmp)
    rSquare <- summary(modeltmp)$r.squared
    dataCorr = rbind(dataCorr, c(nuc,sqrt(rSquare),rSquare,sqrt(mean((subtmp$assigned-subtmp$cs.real)^2)),
                                   mean(abs(subtmp$assigned-subtmp$cs.real))))
  }
  
  
  dataCorr <- dataCorr[rownames(dataCorr) > 1,]
  return(list(h=dataCorr[1,],c=dataCorr[2,],a=dataCorr[3,]))
}