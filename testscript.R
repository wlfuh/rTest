library(pracma)

refdata <- read.table("test.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","error")
subid <- c(1,2,3)

testdata <- refdata[!(refdata$resid %in% subid),"cs"]
moddata <- refdata[(refdata$resid %in% subid),"cs"]
testmatrix <- data.matrix(testdata)
modmatrix <- data.matrix(moddata)
hausdorff_dist(moddata, testmatrix)
stop()

print(hausdorff_dist(testmatrix,jitter(modmatrix, amount = 10)))

# determine the properities of the hausdorff_dist as we data more noise to the subset
# and make plot
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




# try using hungrian algorithm




