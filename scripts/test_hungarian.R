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

# (3) What happens if we repeat this many times? That is, repeat the assignment many time with some noise in moddata, then take the most probable assignment and calculate the cost.
# Is the most probable assignment have an error near zero?
# William: do you want to write the code to test this out?
# How would you accumulate the most assignment?
# How would you determine the most probable?
# Note that the assign function, which you can find in assignments.R, returns both the vector of assignment and an assignment matrix.
# I think the latter would be helpful

# Repeat assignment, with same noise level, get standard deviation of data
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

plot_noise <- function(noiselevel=100){
  data <- NULL
  seq <- NULL
  stdev <- NULL
  for(i in 1:noiselevel){
    tempdata <- getspread(noise=i)
    stdev <- c(stdev, sd(tempdata))
    data <- c(data, tempdata)
    seq <- c(seq, rep(i, length(tempdata)))
  }
  layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
  plot(seq, data, xlab="Noise Level", ylab="Cost", main="Cost Values by Noise Level")
  plot(1:noiselevel, stdev, xlab="Noise Level", ylab="Standard Deviation", main="Standard Deviation by Noise Level (100 iterations)")
  lines(1:noiselevel, stdev, pch=16)
}