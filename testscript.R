library(pracma)

refdata <- read.table("test.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","error")
subid <- c(1,2,3)

plot_haus <- function(iter = 100){
  testdata <- refdata[!(refdata$resid %in% subid),"cs"]
  moddata <- refdata[(refdata$resid %in% subid),"cs"]
  testmatrix <- data.matrix(testdata)
  modmatrix <- data.matrix(moddata)
  hausdorff_dist(moddata, testmatrix)
  
  print(hausdorff_dist(testmatrix,jitter(modmatrix, amount = 10)))
  
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
  
  plot(ms, final_dist, main="Hausdorff Distance depending on Jitter Amount", 
       xlab="Jitter Amount", ylab="Hausdorff Distance")
}

# try using hungrian algorithm

plot_hung <- function(iter = 100){
  testdata <- refdata[!(refdata$resid %in% subid),"cs"]
  moddata <- refdata[(refdata$resid %in% subid),"cs"]
  print(find_cost(moddata, testdata))
  final_dist <- NULL
  ms <- seq(1,iter,1)
  for (m in ms){
    # average distance
    tmp_dist <- NULL
    for (i in 1:iter){
      tmp_dist <- c(tmp_dist, find_cost(testdata,jitter(moddata, amount = m)))
    }
    final_dist <- c(final_dist, mean(tmp_dist))
  }
  plot(ms, final_dist, main="Hungarian Algorithm Cost depending on Jitter Amount", 
       xlab="Jitter Amount", ylab="Total Cost")
}

plot_special <- function(iter = 100){
  testdata <- refdata[!(refdata$resid %in% subid),"cs"]
  moddata <- refdata[(refdata$resid %in% subid),"cs"]
  print(find_cost(moddata, testdata))
  final_dist <- NULL
  ms <- seq(1,iter,1)
  for (m in ms){
    # average distance
    tmp_dist <- NULL
    for (i in 1:iter){
      modjitt <- jitter(moddata, amount = m)
      tmp_dist <- c(tmp_dist, mean(find_cost(testdata,modjitt), hausdorff_dist(testdata,modjitt)))
    }
    final_dist <- c(final_dist, mean(tmp_dist))
  }
  plot(ms, final_dist, main="Average cost/distance (Hungarian and Hausdorff combined)", 
       xlab="Jitter Amount", ylab="Cost/ Distance")
  
}

find_cost <- function(x, y){
  require(clue)
  # maps chemical shifts in y to chemical shifts in x 
  # Input -- x (vector): chemical shifts 
  # Input -- y (vector): chemical shifts which are to be mapped to those in x (not x can be greater than y)
  if(length(x) < length(y)){
    temp = x
    x = y
    y = temp
  }
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
  costmat <- abs(xmat-ymat)
  a <- solve_LSAP(costmat)
  cost = as.integer(0)
  for(i in 1:length(a))
    cost = cost + costmat[i, a[i]]
  return(cost)
  # cbind(x[assignments],y)
}

# k- nearest algorithm
# uses library class
plot_knn <- function(noise = 0){
  standData <- read.table("standard.dat", header=FALSE)
  colnames(standData) <- c("nuc", "atom", "aType", "cs", "minCs", "maxCs", "avgCs", "sd")
  train <- data.matrix(refdata[!(refdata$resid %in% subid),"cs"])
  test <- data.matrix(refdata[(refdata$resid %in% subid),"cs"])
  test <- jitter(test, amount=noise)
  # classification by residue name
  cl <- factor(refdata[!(refdata$resid %in% subid),"resname"])
  kvals <- seq(1,length(train),1)
  diffs <- NULL
  # unsure if proper classification/ way to find different
  for(i in 1:length(train)){
    knnres <- knn(train, test, cl, k=i)
    knnmatrix <- data.matrix(knnres)
    originClass <- data.matrix(refdata[(refdata$resid %in% subid),"resname"])
    numDiffs <- 0
    if(length(knnmatrix) != length(originClass))
      stop("KNN result has different length than the original classification of test")
    for(j in 1:length(knnres)){
      if(knnmatrix[j] != originClass[j])
        numDiffs <- numDiffs + 1
    }
    diffs <- c(diffs, numDiffs)
  }
  plot(kvals, diffs)
}

# linear regression
plot_linear <- function(){
  testdata <- refdata[!(refdata$resid %in% subid),"cs"]
  moddata <- refdata[(refdata$resid %in% subid),"cs"]
  test <- data.frame(id=1:length(testdata),cs=testdata)
  mod <- data.frame(id=1:length(moddata),cs=moddata)
  #plot(seq(1,length(testmatrix),1), testmatrix)
  #fit <- lm(id ~ cs, data=test)
  #plot(fit)
  lm.out = lm(cs ~ id, data=test)
  lm.out2 = lm(cs ~ id, data=mod)
  layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
  plot(cs ~ id, data=test, main="Test Plot")
  abline(lm.out, col="red")
  plot(cs ~ id, data=mod, main="Mod Plot")
  abline(lm.out2, col="green")
  print(coefficients(lm.out))
  print(coefficients(lm.out2))
}

'
knnres <- knn(train, test, cl, k=7)
print(knnres)
'
'
knnmatrix <- data.matrix(knnres)
originClass <- data.matrix(refdata[(refdata$resid %in% subid),"resname"])
numDiffs <- 0
if(length(knnmatrix) != length(originClass))
  stop("KNN result has different length then original classification of test")
for(j in 1:length(knnres)){
  if(knnmatrix[j] == originClass[j])
    numDiffs <- numDiffs + 1
}
'

'
testKn = knn(modmatrix, testmatrix, modmatrix, k=30)
input_testKn <- data.frame(testKn)
plot(input_testKn)
'

'
train <- rbind(iris3[1:25,,1], iris3[1:25,,2], iris3[1:25,,3])
test <- rbind(iris3[26:50,,1], iris3[26:50,,2], iris3[26:50,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
knn(train, test, cl, k = 3, prob=TRUE)
attributes(.Last.value)
'


