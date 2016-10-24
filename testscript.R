testdata <- read.table("test.dat", header=FALSE)
moddata <- read.table("test_modified.dat", header=FALSE)

testraw <- testdata[4]

modraw <- moddata[4]

testmatrix <- data.matrix(testraw)
modmatrix <- data.matrix(modraw)
print(hausdorff_dist(testmatrix, modmatrix))
# install.packages : pracma