# goto working directory
setwd("~/GitSoftware/rTest/")

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

