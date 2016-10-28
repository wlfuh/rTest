library(pracma)

refdata <- read.table("test.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","error")
subid <- c(1,2,3)

testdata <- refdata[!(refdata$resid %in% subid),"cs"]
moddata <- refdata[(refdata$resid %in% subid),"cs"]
testmatrix <- data.matrix(testdata)
modmatrix <- data.matrix(moddata)
hausdorff_dist(moddata, testmatrix)

print(hausdorff_dist(testmatrix,jitter(modmatrix, amount = 10)))

# determine the properities of the hausdorff_dist as we data more noise to the subset
# and make plot
"
final_dist <- NULL
ms <- seq(1,100,1)
for (m in ms){
  # average distance
  tmp_dist <- NULL
  for (i in 1:100){
    tmp_dist <- c(tmp_dist, hausdorff_dist(testmatrix,jitter(modmatrix, amount = m)))
  }
  final_dist <- c(final_dist, mean(tmp_dist))
}

plot(ms, final_dist)
"
# try using hungrian algorithm

# k- nearest algorithm
testKn = knn(modmatrix, testmatrix, modmatrix, k=30)
input_testKn <- data.frame(testKn)
plot(input_testKn)
'
train <- rbind(iris3[1:25,,1], iris3[1:25,,2], iris3[1:25,,3])
test <- rbind(iris3[26:50,,1], iris3[26:50,,2], iris3[26:50,,3])
cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
knn(train, test, cl, k = 3, prob=TRUE)
attributes(.Last.value)
'


