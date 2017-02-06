rmsd <- read.table("data/rmsd.txt", col.names = c("state","rmsd"))

# helper function for analyze_x
getModelError <- function(modeldata, type){
  if(type == "actual")
    cmpdata <- modeldata$actual
  else if(type == "predicted")
    cmpdata <- modeldata$cs
  else
    stop("Invalid type")
  
  errors <- data.frame(
    state = modeldata$state[1],
    mae = mean(abs(modeldata$assigned-cmpdata)),
    rmse = sqrt(mean((modeldata$assigned-cmpdata)^2)),
    r = cor(modeldata$assigned, cmpdata, method = "pearson"),
    rho = cor(modeldata$assigned, cmpdata, method = "spearman"),
    tau = cor(modeldata$assigned, cmpdata, method = "kendall")
  )
}


# assigned vs actual cs
analyze_native <- function(filename){
  require(plyr)
  modeldata <- read.table(filename, header=TRUE)
  
  modelerror <- ddply(.dat=modeldata, .var=c("state"), 
                      .fun = getModelError, type="actual")
  
  return(merge(modelerror, rmsd))
}


# assigned vs predicted cs
analyze_pred <- function(filename){
  require(plyr)
  modeldata <- read.table(filename, header=TRUE)
  
  modelerror <- ddply(.dat=modeldata, .var=c("state"), 
                      .fun = getModelError, type="predicted")
  
  return(merge(modelerror, rmsd))
}

# return rmsd error for each error method
get_rmsd_error <- function(result){
  colsRes <- colnames(result)  
  
  out <- data.frame( mae=0, rmse=0, r=0, rho=0,tau=0, rmsd = 0)
  
  trueOrder <- rmsd[order(rmsd$rmsd),]
  
  for(col in colsRes){
    if(col == 'state')
      next
    if(col == 'r' || col == 'rho' || col == 'tau')
      out[col] <- cor(result[order(result[col], decreasing = TRUE),]$rmsd, trueOrder$rmsd, method="spearman")    
    else
      out[col] <- cor(result[order(result[col]),]$rmsd, trueOrder$rmsd, method="spearman")  
  }
  return(out)
}
