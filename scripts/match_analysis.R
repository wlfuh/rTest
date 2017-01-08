# load library and our custom functions
library(pracma)
library(plyr)
source("scripts/hung_trials.R")

# read in reference data
refdata <- read.table("data/measured_shifts_2KOC.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","dummy")

# model predicted chemical shift
modeldata <- read.table("data/predicted_CS_table_test_clean.txt", header=TRUE)

# larmord accuracy data
accu <- read.table("data/larmord_accuracy_resname_nucleus.txt",header=FALSE)
colnames(accu) <- c("resname", "nucleus", "error")

matchModel <- function(){
  levels <- length(unique(modeldata$state))
  
  for(i in 1:levels){
    predcs <- subset(modeldata, modeldata$state==i)
    predcs <- rename(predcs, c("nuclei" = "nucleus", "predCS"="cs"))
    predcs <- merge(predcs, accu)
    predcs <- predcs[,c(1,4,2,5,6,3)]
    predcs <- predcs[order(predcs$resid, predcs$nucleus),]
    print(hung_trials(predcs))
  }
  
}