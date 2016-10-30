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
  
  plot(ms, final_dist)
}

# try using hungrian algorithm

# k- nearest algorithm
# uses library class
plot_knn <- function(noise = 0){
  standData <- read.table("standard.dat", header=FALSE)
  colnames(standData) <- c("nuc", "atom", "aType", "cs", "minCs", "maxCs", "avgCs", "sd")
  train <- data.matrix(refdata[!(refdata$resid %in% subid),"cs"])
  test <- data.matrix(refdata[(refdata$resid %in% subid),"cs"])
  test <- jitter(test, amount=noise)
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
      stop("KNN result has different length then original classification of test")
    for(j in 1:length(knnres)){
      if(knnmatrix[j] != originClass[j])
        numDiffs <- numDiffs + 1
    }
    diffs <- c(diffs, numDiffs)
  }
  plot(kvals, diffs)
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


