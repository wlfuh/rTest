# load library and our custom functions
library(pracma)
library(plyr)
source("scripts/assignments.R")

# read in reference data
refdata <- read.table("data/measured_shifts_2KOC.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","dummy")

refdata2 <- read.table("data/test.dat", header=FALSE)
colnames(refdata2) <- c("resname","resid","nucleus","cs","error")

# vector of reference values
testdata <- refdata[,"cs"]
testdata2 <- refdata2[,"cs"]

# vector of subset
predcs <- read.table("data/predicted_shifts_2KOC.dat", header=FALSE)
colnames(predcs) <- c("resname","resid","nucleus","cs","dummy")
accu <- read.table("data/larmord_accuracy_resname_nucleus.txt",col.names = c("resname","nucleus","error"))
predcs <- merge(predcs, accu)


moddata <- predcs$cs


# (1) Sanity check: given that moddata is simply a copy data contained in testdata, then assignment error should be 0
a <- assign(testdata,moddata)
print(get_assignment_cost(a$a, a$costmat))

# (2) What happen if we add some random noise?
a <- assign(testdata,jitter(moddata))
print(get_assignment_cost(a$a, a$costmat))

# Numerical Summary
print(summary(testdata))
print("Standard Deviation")
print(sd(testdata))

add_noise <- function(x, scale=1){
  x$cs+rnorm(1, 0, scale*x$error)
}

getspread <- function(iter=100, noise=1, diff=NULL){
  data <- NULL
  for(i in 1:iter){
    a <- assign(testdata,jitter(moddata, factor=noise, amount=diff))
    data <- c(data, get_assignment_cost(a$a, a$costmat))
  }
  # print(summary(data))
  #print(paste("standard deviation", sd(data)))
  return(data)
}

# ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=0.1)

# Plots Assignment Cost based on scaled noise level, up to maxscale, split into inc times each iter times
# comp- TRUE- plot w/ both refdata, FALSE- plot w/ only one refdata
plot_scaled <- function(maxscale=1, inc, iter=1, doplot=TRUE){
  costs <- NULL
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:inc)
  filename <- paste("output/plot_scaled_output_", inc, "_increment_", iter, "_times.txt")
  bestassign <- assign(testdata, predcs$cs)
  bestcost <- get_assignment_cost(bestassign$a, bestassign$costmat)
  bestscale <- 0
  # stdevs <- NULL
  for(i in scales){
    ith_cost <- NULL
    for(j in 1:iter){
      predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
      a <- assign(testdata, predcs_mod$V1)
      pred_cost = get_assignment_cost(a$a, a$costmat)
      if(pred_cost < bestcost){
        bestassign <- a
        bestcost <- pred_cost
        bestscale <- i
      }
      ith_cost <- c(ith_cost, as.numeric(pred_cost))
    }
    costs <- c(costs, mean(ith_cost))
    #if(iter > 1)
      #stdevs <- c(stdevs,sd(ith_cost))
  }
  if(doplot)
    plot(scales, costs, xlab=paste("Scale (incrementing by ", maxscale/inc, " )"), ylab="Hungarian Assignment Cost",
       main=paste("Hungarian Assignment Costs based on Scaled Random Norm Errors\n", iter, " iterations for each"))
  print(summary(costs))
  print("Best")
  print(bestassign$a)
  print(paste("Cost: ",bestcost))
  print(paste("Scale: ",bestscale))
}

# Things to test, different maxscale, different iterations
# Multiple Iterations at each scale level

# calculates number of differences b/w two LSAP (assumed to be same size)
get_diffs <- function(l1, l2){
  diff = 0
  for(i in 1:length(l1)){
    if(l1[i] != l2[i])
      diff <- diff + 1
  }
  return(diff)
}

# plot the number of different assignment at each inc relative to original assignments
plot_diffs <- function(maxscale=1, inc){
  original <- assign(testdata, predcs$cs)
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:inc)
  diffs <- 0
  for(i in scales){
    if(i == 0)
      next
    predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
    a <- assign(testdata, predcs_mod$V1)
    diffs <- c(diffs, get_diffs(original$a, a$a))
  }
  plot(scales, diffs, xlab=paste("Scale (incrementing by ", maxscale/inc, " )"), ylab="Differences",
       main="Number of Different Assignment at each scale")
}

hung_trials <- function(maxscale, inc, iter){
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:inc)
  masterlist <- list()
  for(i in scales){
    if(i == 0)
      next
    mat <- matrix(NA, nrow=iter, ncol=nrow(predcs))
    prob <- matrix(0, nrow=nrow(predcs), ncol=nrow(refdata))
    scale_costs <- NULL
    acc_resid <- NULL
    acc_nucleus <- NULL
    
    for(j in 1:iter){
      predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
      a <- assign(testdata, predcs_mod$V1)
      pred_cost = get_assignment_cost(a$a, a$costmat)
      scale_costs <- c(scale_costs, pred_cost)
      
      mat[j,] <- as.vector(a$a)
      prob = prob + a$a_mat
      d <- cbind(predcs_mod, refdata[as.vector(a$a),])
      acc_resid <- c(acc_resid, mean(as.character(d[,2])==as.character(d[,6])))
      acc_nucleus <- c(acc_nucleus, mean(d[,1]==d[,5]))
    }
    key <- paste("scale",i,sep="")
    masterlist[[key]] = list(costs=scale_costs, assignments=mat, iterations=iter, probs=prob/iter, 
                             acc_resid=acc_resid, acc_nuc=acc_nucleus)
  }
  save(masterlist, file=paste("output/",maxscale,"scale_",(maxscale/(inc - 1)),"inc_",iter,"iter.RData",sep=""))
}

