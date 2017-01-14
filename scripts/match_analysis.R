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

matchModel <- function(iter=5000, scale=0.75,method="mean"){
  levels <- length(unique(modeldata$state))
  
  hData <- data.frame(model=0,nuc="H",r=0,rsquare=0,RMSE=0,MAE=0,stringsAsFactors=FALSE)
  cData <- data.frame(model=0,nuc="C",r=0,rsquare=0,RMSE=0,MAE=0,stringsAsFactors=FALSE)
  #print(paste("output/hTest",iter,"_iter_",scale,"_scale.csv",sep=""))
  for(i in 1:levels){
    predcs <- subset(modeldata, modeldata$state==i)
    predcs <- rename(predcs, c("nuclei" = "nucleus", "predCS"="cs"))
    predcs <- merge(predcs, accu)
    predcs <- predcs[,c(1,4,2,5,6,3)]
    predcs <- predcs[order(predcs$resid, predcs$nucleus),]
    
    # TODO save results in larger table one for H and C, save table to csv
    data <- hung_trials(predcs, iter, scale)
    
    data$h["model"] <- c(i)
    data$c["model"] <- c(i)
    
    data$h <- data$h[,c(6,1,2,3,4,5)]
    data$c <- data$c[,c(6,1,2,3,4,5)]
    
    hData = rbind(hData, data$h)
    cData = rbind(cData, data$c)
  }
  
  write.csv(hData[rownames(hData) > 1,], 
            file=paste("output/hTest_",iter,"_iter_",scale,"_scale_",method,"_method.csv",sep=""), row.names = FALSE)
            #, col.names = c("model","nuc","r","rsquare","RMSE","MAE"))
  write.csv(cData[rownames(cData) > 1, ], 
            file=paste("output/cTest",iter,"_iter_",scale,"_scale_",method,"_method.csv",sep=""), row.names = FALSE)
            #, col.names = c("model","nuc","r","rsquare","RMSE","MAE"))
  
}