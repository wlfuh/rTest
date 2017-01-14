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
  predcs$assigned <- refdata[a,"cs"]
  
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
get_correlation <- function(res, nuclei=c("H","C")){
  
  tmp <- merge(res$cs,refdata, by=c("resid","nucleus"),suffixes=c(".pred",".real"))
  dataCorr <- data.frame(nuc="B1",r=0,rsquare=0,RMSE=0,MAE=0,stringsAsFactors=FALSE)
  
  for(nuc in nuclei){
    subtmp <- tmp[substr(tmp$nucleus,1,1)==nuc,]
    modeltmp <- lm(assigned~cs.real, data=subtmp)
    rSquare <- summary(modeltmp)$r.squared
    dataCorr = rbind(dataCorr, c(nuc,sqrt(rSquare),rSquare,sqrt(mean((subtmp$assigned-subtmp$cs.real)^2)),
                                   mean(abs(subtmp$assigned-subtmp$cs.real))))
  }
  
  dataCorr <- dataCorr[rownames(dataCorr) > 1,]
  return(list(h=dataCorr[1,],c=dataCorr[2,]))
}