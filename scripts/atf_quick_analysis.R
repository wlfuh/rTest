analyze_1 <- function(filename){
  require(plyr)
  modeldata <- read.table(filename, header=TRUE)
  
  modelerror <- ddply(.dat=modeldata, .var=c("state"), 
                      .fun = function(modeldata){
                        errors <- data.frame(state=0, mae=0, rmse=0, r=0)
                        errors$state <- modeldata$state[1]
                        # TODO MAE, RMSE
                        errors$mae <-  mean(abs(modeldata$assigned-modeldata$actual))
                        errors$rmse <- sqrt(mean((modeldata$assigned-modeldata$actual)^2))
                        errors$r <- cor(modeldata$assigned, modeldata$actual, method = "pearson")
                        errors$rho <- cor(modeldata$assigned, modeldata$actual, method = "spearman")
                        errors$tau <- cor(modeldata$assigned, modeldata$actual, method = "kendall")
                        errors
                      })
  return(modelerror)
}


analyze_2 <- function(filename){
  require(plyr)
  modeldata <- read.table(filename, header=TRUE)
  
  modelerror <- ddply(.dat=modeldata, .var=c("state"), 
                      .fun = function(modeldata){
                        errors <- data.frame(state=0, mae=0, rmse=0, r=0)
                        errors$state <- modeldata$state[1]
                        # TODO MAE, RMSE
                        errors$mae <-  mean(abs(modeldata$assigned-modeldata$cs))
                        errors$rmse <- sqrt(mean((modeldata$assigned-modeldata$cs)^2))
                        errors$r <- cor(modeldata$assigned, modeldata$cs, method = "pearson")
                        errors$rho <- cor(modeldata$assigned, modeldata$cs, method = "spearman")
                        errors$tau <- cor(modeldata$assigned, modeldata$cs, method = "kendall")
                        errors
                      })
  return(modelerror)
}


# Q1: Are the models that exhibit the lowest assignment error native-like, i.e., have low RMSD values? 
# This is important for validating structure-based assignment.
rest <- analyze_1("output/qm/predicted_shifts_qm_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$mae),]

rest <- analyze_1("output/ramsey/predicted_shifts_ramsey_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$mae),]

rest <- analyze_1("output/larmord/predicted_shifts_larmord_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$mae),]

rest <- analyze_1("output/combined/predicted_shifts_combined_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$rmse),]


# Q2: Can we the error between assigned chemical shift and predicted chemical shift be used to identify the native structure? 
# This is important for demonstrating its utility
rest <- analyze_2("output/qm/predicted_shifts_qm_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$mae),]

rest <- analyze_2("output/ramsey/predicted_shifts_ramsey_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$mae),]

rest <- analyze_2("output/larmord/predicted_shifts_larmord_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$mae),]

rest <- analyze_2("output/combined/predicted_shifts_combined_modified.txt")
rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))
rest <- merge(rest, rmsd)
rest[order(rest$rmse),]





