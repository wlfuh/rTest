# load library and our custom functions
library(pracma)
source("scripts/assignments.R")

# read in reference data
refdata <- read.table("data/test.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","error")
subid <- c(1,2,3)

# vector of reference values
testdata <- refdata[,"cs"]
# vector of subset
moddata <- refdata[(refdata$resid %in% subid),"cs"]

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

# (3) What happens if we repeat this many times? That is, repeat the assignment many time with some noise in moddata, then take the most probable assignment and calculate the cost.
# Is the most probable assignment have an error near zero?
# William: do you want to write the code to test this out?
# How would you accumulate the most assignment?
# How would you determine the most probable?
# Note that the assign function, which you can find in assignments.R, returns both the vector of assignment and an assignment matrix.
# I think the latter would be helpful

# Need to determine the maximum cost, to be considered probable assignment, perhaps by standard deviation of entire data
# But the range of values differs depending on the nuclues
# For instance, H ranges from 3.856 to 14.458 and jumps to 96.6 for C
# refdata[with(refdata, order(cs)), ]

# Could also use mean, but does not really make sense.


# Gets assignment costs for i iterations, at noise level n
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

# Plots Assignment Cost as a function of noiselevel, with each noiselevel plotted i iter times
plot_noise <- function(noiselevel=100, iter=100, table=FALSE){
  data <- NULL
  seq <- NULL
  stdev <- NULL
  pvals <- NULL
  for(i in 1:noiselevel){
    tempdata <- getspread(iter,noise=i)
    stdev <- c(stdev, sd(tempdata))
    pvals <- c(pvals, t.test(tempdata)$p.value)
    data <- c(data, tempdata)
    seq <- c(seq, rep(i, length(tempdata)))
  }
  if(!table){
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
    plot(seq, data, xlab="Noise Level", ylab="Cost", main="Cost Values by Noise Level")
    plot(1:noiselevel, stdev, xlab="Noise Level", ylab="Standard Deviation", main=paste("Standard Deviation by Noise Level (",iter,"iterations)",sep=" "))
    lines(1:noiselevel, stdev, pch=16)
    # plot(1:noiselevel, pvals, xlab="Noise Level", ylab="p-value", main=paste("p-value by Noise Level (",iter,"iterations)",sep=" "))
    # lines(1:noiselevel, pvals, pch=16)
  }
  else{
    print(cbind(seq, data))
  }
}

# Plots Assignment Cost as a function of noiselevel, averaging i iterations of the assignment cost 
# at each noise level
plot_noise_avg <- function(noiselevel=100, iter=100, table=FALSE){
  data <- NULL
  stdev <- NULL
  # pvals <- NULL
  for(i in 1:noiselevel){
    tempdata <- getspread(iter,noise=i)
    stdev <- c(stdev, sd(tempdata))
    # pvals <- c(pvals, t.test(tempdata)$p.value)
    data <- c(data, mean(tempdata))
  }
  if(!table){
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
    plot(1:noiselevel, data, xlab="Noise Level", ylab="Cost", main="Cost Values by Noise Level (Averaging Iterations)")
    dataMean = mean(testdata)
    abline(h=dataMean)
    text(noiselevel, dataMean, paste("mean ",dataMean), pos="3")
    plot(1:noiselevel, stdev, xlab="Noise Level", ylab="Standard Deviation", main=paste("Standard Deviation by Noise Level (",iter,"iterations)",sep=" "))
    lines(1:noiselevel, stdev, pch=16)
    # plot(1:noiselevel, pvals, xlab="Noise Level", ylab="p-value", main=paste("p-value by Noise Level (",iter,"iterations)",sep=" "))
    # lines(1:noiselevel, pvals, pch=16)
  }
  else{
    print(cbind(1:noiselevel, data))
  }
}

# Plot assingment cost as a function of the number iteration at a fixed noise level
plot_iter <- function(iter=100){
  data <- NULL
  seq <- NULL
  stdev <- NULL
  for(i in 1:iter){
    tempdata <- getspread(iter=i)
    stdev <- c(stdev, sd(tempdata))
    data <- c(data, tempdata)
    seq <- c(seq, rep(i, length(tempdata)))
  }
  layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
  plot(seq, data, xlab="Iterations", ylab="Cost", main="Cost Values by Iterations")
  plot(1:iter, stdev, xlab="Iterations", ylab="Standard Deviation", main="Standard Deviation by Iterations (Noise Level 1)")
  lines(1:iter, stdev, pch=16)
}

test_reverse <- function(){
  revmod <- rev(moddata)
  a <- assign(testdata,revmod)
  print(get_assignment_cost(a$a, a$costmat))
  
}

test_sample <- function(){
  sampmod <- sample(moddata)
  a <- assign(testdata,sampmod)
  print(get_assignment_cost(a$a, a$costmat))
  
}

