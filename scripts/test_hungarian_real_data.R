# load library and our custom functions
library(pracma)
library(plyr)
source("scripts/assignments.R")

# read in reference data
refdata <- read.table("data/measured_shifts_2KOC.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","dummy")

# vector of reference values
testdata <- refdata[,"cs"]

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
plot_scaled <- function(maxscale=1, inc, iter=1){
  costs <- NULL
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:inc)
  stdevs <- NULL
  for(i in scales){
    ith_cost <- NULL
    for(j in 1:iter){
      predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
      a <- assign(testdata, predcs_mod$V1)
      pred_cost = get_assignment_cost(a$a, a$costmat)
      ith_cost <- c(ith_cost, as.numeric(pred_cost))
    }
    costs <- c(costs, median(ith_cost))
    if(iter > 1)
      stdevs <- c(stdevs,sd(ith_cost))
  }
  plot(scales, costs, xlab=paste("Scale (incrementing by ", maxscale/inc, " )"), ylab="Hungarian Assignment Cost",
       main=paste("Hungarian Assignment Costs based on Scaled Random Norm Errors\n", iter, " iterations for each"))
  if(iter > 1){
    stdframe <- data.frame(scales, stdevs)
    names(stdframe) <- c("Scale", "Standard Deviation")
    print(stdframe)
  }
}

# Things to test, different maxscale, different iterations
# Multiple Iterations at each scale level


