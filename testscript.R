library(pracma)

testdata <- read.table("test.dat", header=FALSE)
moddata <- read.table("test_modified.dat", header=FALSE)

testraw <- testdata[4]

modraw <- moddata[4]
testmatrix <- data.matrix(testraw)
modmatrix <- data.matrix(modraw)
print(hausdorff_dist(testmatrix,jitter(modmatrix, amount = 10)))
# install.packages : pracma


# determine the properities of the hausdorff_dist as we data more noise to the subset
final_dist <- NULL
for (m in seq(1,100,1)){
  # average distance
  tmp_dist <- NULL
  for (i in 1:100){
    tmp_dist <- c(tmp_dist, hausdorff_dist(testmatrix,jitter(modmatrix, amount = m)))
  }
}


