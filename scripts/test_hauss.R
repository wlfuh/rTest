library(pracma)

refdata <- read.table("data/test.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","error")
subid <- c(1,2,3)

# vector of reference values
testdata <- refdata[,"cs"]
# vector of subset
moddata <- refdata[(refdata$resid %in% subid),"cs"]

testmatrix <- data.matrix(testdata)
modmatrix <- data.matrix(moddata)
print(hausdorff_dist(moddata, testmatrix))

plot_haus <- function(iter = 100, axises=TRUE){

  
  #print(hausdorff_dist(testmatrix,jitter(modmatrix, amount = 10)))
  
  # determine the properities of the hausdorff_dist as we data more noise to the subset
  # and make plot
  
  final_dist <- NULL
  ms <- seq(1,iter,1)
  for (m in ms){
    # average distance
    tmp_dist <- NULL
    for (i in 1:iter){
      tmp_dist <- c(tmp_dist, hausdorff_dist(testmatrix,jitter(modmatrix, amount = m)))
    }
    final_dist <- c(final_dist, mean(tmp_dist))
  }
  if(axises)
    plot(ms, final_dist, main="Hausdorff Distance depending on Jitter Amount", 
         xlab="Jitter Amount", ylab="Hausdorff Distance")
  else
    points(ms, final_dist, xlab="", ylab="")
  print(capture.output(cat("Hausdorff Standard Deviation: ", sd(final_dist))))
  print(capture.output(cat("Hausdorff Mean: ", mean(final_dist))))
}

test_reverse <- function(){
  revmod <- rev(moddata)
  print(hausdorff_dist(revmod, moddata))
  print(hausdorff_dist(moddata, testmatrix))
  print(hausdorff_dist(revmod, testmatrix))
}

test_sample <- function(){
  sampmod <- sample(moddata)
  print(hausdorff_dist(sampmod, moddata))
  print(hausdorff_dist(moddata, testmatrix))
  print(hausdorff_dist(sampmod, testmatrix))
}

# Haussdorff distance can be used to see how much a subset matches with original set

plot_sub <- function(noise=0){
  x <- as.data.frame(table(refdata[,"resid"]))
  resids <- 1:length(x[,"Var1"])
  data <- NULL
  for(i in resids){
    submod <- jitter(refdata[(refdata$resid %in% 1:i),"cs"], factor=noise)
    #print(refdata[(refdata$resid %in% 1:i),"cs"])
    #print(submod)
    data <- c(data, hausdorff_dist(data.matrix(submod), testmatrix))
  }
  plot(resids, data)
  print(data)
}