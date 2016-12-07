# load library and our custom functions
library(pracma)
library(plyr)
source("scripts/assignments.R")

# read in reference data
refdata <- read.table("data/measured_shifts_2KOC.dat", header=FALSE)
colnames(refdata) <- c("resname","resid","nucleus","cs","dummy")

refdata2 <- read.table("data/test.dat", header=FALSE)
colnames(refdata2) <- c("resname","resid","nucleus","cs","error")

# vector of reference values
testdata <- refdata[,"cs"]
testdata2 <- refdata2[,"cs"]

# vector of subset
predcs <- read.table("data/predicted_shifts_2KOC.dat", header=FALSE)
colnames(predcs) <- c("resname","resid","nucleus","cs","dummy")
accu <- read.table("data/larmord_accuracy_resname_nucleus.txt",col.names = c("resname","nucleus","error"))
predcs <- merge(predcs, accu)
predcs <- predcs[order(predcs$resid, predcs$nucleus),]


moddata <- predcs$cs


# (1) Sanity check: given that moddata is simply a copy data contained in testdata, then assignment error should be 0
a <- assign(testdata,moddata)
print(get_assignment_cost(a$a, a$costmat))

# (2) What happen if we add some random noise?
a <- assign(testdata,jitter(moddata))
print(get_assignment_cost(a$a, a$costmat))

# Numerical Summary
print(summary(testdata))
print("Standard Deviation")
print(sd(testdata))

add_noise <- function(x, scale=1){
  x$cs+rnorm(1, 0, scale*x$error)
}

getspread <- function(iter=100, noise=1, diff=NULL){
  data <- NULL
  for(i in 1:iter){
    a <- assign(testdata,jitter(moddata, factor=noise, amount=diff))
    data <- c(data, get_assignment_cost(a$a, a$costmat))
  }
  # print(summary(data))
  #print(paste("standard deviation", sd(data)))
  return(data)
}

# ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=0.1)

# Plots Assignment Cost based on scaled noise level, up to maxscale, split into inc times each iter times
# comp- TRUE- plot w/ both refdata, FALSE- plot w/ only one refdata
plot_scaled <- function(maxscale=1, inc, iter=1, doplot=TRUE){
  costs <- NULL
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:inc)
  filename <- paste("output/plot_scaled_output_", inc, "_increment_", iter, "_times.txt")
  bestassign <- assign(testdata, predcs$cs)
  bestcost <- get_assignment_cost(bestassign$a, bestassign$costmat)
  bestscale <- 0
  sink(filename)
  # stdevs <- NULL
  for(i in scales){
    ith_cost <- NULL
    for(j in 1:iter){
      predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
      predcs_mod <- predcs_mod[order(predcs_mod$resid, predcs_mod$nucleus),]
      a <- assign(testdata, predcs_mod$V1)
      pred_cost = get_assignment_cost(a$a, a$costmat)
      if(pred_cost < bestcost){
        bestassign <- a
        bestcost <- pred_cost
        bestscale <- i
      }
      ith_cost <- c(ith_cost, as.numeric(pred_cost))
      print(a$a)
    }
    costs <- c(costs, mean(ith_cost))
    #if(iter > 1)
    #stdevs <- c(stdevs,sd(ith_cost))
  }
  if(doplot)
    plot(scales, costs, xlab=paste("Scale (incrementing by ", maxscale/inc, " )"), ylab="Hungarian Assignment Cost",
         main=paste("Hungarian Assignment Costs based on Scaled Random Norm Errors\n", iter, " iterations for each"))
  sink(type = "message")
  sink()
  print(summary(costs))
  print("Best")
  print(bestassign$a)
  print(paste("Cost: ",bestcost))
  print(paste("Scale: ",bestscale))
}

# Perform multiple iterations at each scale and stores an R object with data in it
hung_trials <- function(maxscale, inc, iter, minscale=0,scales=NULL, notes="",savef=TRUE){
  
  if(inc == 1 && is.null(scales))
    scales = maxscale
  else{
    range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale + minscale }
    scales <- range01(1:inc)  
  }
  
  masterlist <- list()
  for(i in scales){
    if(i == 0)
      next
    mat <- matrix(NA, nrow=iter, ncol=nrow(predcs))
    prob <- matrix(0, nrow=nrow(predcs), ncol=nrow(refdata))
    scale_costs <- NULL
    acc_resid <- NULL
    acc_nucleus <- NULL
    
    for(j in 1:iter){
      predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
      predcs_mod <- predcs_mod[order(predcs_mod$resid, predcs_mod$nucleus),]
      a <- assign(testdata, predcs_mod$V1)
      pred_cost = get_assignment_cost(a$a, a$costmat)
      scale_costs <- c(scale_costs, pred_cost)
      
      mat[j,] <- as.vector(a$a)
      prob = prob + a$a_mat
      d <- cbind(predcs_mod, refdata[as.vector(a$a),])
      acc_resid <- c(acc_resid, mean(as.character(d[,2])==as.character(d[,6])))
      acc_nucleus <- c(acc_nucleus, mean(d[,1]==d[,5]))
    }
    key <- paste("scale",i,sep="")
    masterlist[[key]] = list(costs=scale_costs, assignments=mat, iterations=iter, probs=prob/iter, 
                             acc_resid=acc_resid, acc_nuc=acc_nucleus)
  }
  if(savef){
    if(inc == 1)
      scalep <- maxscale
    else
      scalep <- (maxscale/(inc - 1))
    filename <- paste("output/",notes,maxscale,"scale_",scalep,"inc_",iter,"iter.RData",sep="")
    save(masterlist, file=filename)
  }
}


# Things to test, different maxscale, different iterations
# Multiple Iterations at each scale level

# calculates number of differences b/w two LSAP (assumed to be same size)
get_diffs <- function(l1, l2){
  diff = 0
  for(i in 1:length(l1)){
    if(l1[i] != l2[i])
      diff <- diff + 1
  }
  return(diff)
}

# plot the number of different assignment at each inc relative to original assignments
plot_diffs <- function(maxscale=1, inc){
  original <- assign(testdata, predcs$cs)
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:inc)
  diffs <- 0
  for(i in scales){
    if(i == 0)
      next
    predcs_mod <- ddply(.dat=predcs, .var=c("resid","nucleus"), .fun = add_noise, scale=i)
    
    a <- assign(testdata, predcs_mod$V1)
    diffs <- c(diffs, get_diffs(original$a, a$a))
  }
  plot(scales, diffs, xlab=paste("Scale (incrementing by ", maxscale/inc, " )"), ylab="Differences",
       main="Number of Different Assignment at each scale")
}

# plot mean accuracy over iterations
plot_accuracy <- function(){
  load("~/rTest/output/1.5scale_0.25inc_10000iter.RData")
  accu_res_total <- NULL
  accu_nuc_total <- NULL
  maxscale <- 1.5
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:(length(masterlist)+1))
  for(i in masterlist){
    accu_res <- i$acc_resid
    accu_nuc <- i$acc_nuc
    accu_res_total <- c(accu_res_total,mean(accu_res))
    accu_nuc_total <- c(accu_nuc_total,mean(accu_nuc))
  }
  layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
  plot(scales[2:length(scales)], accu_res_total, xlab="Scale", ylab="Accuracy", main=paste("Accuracy of Residue ID by Scale",
                                                                                           masterlist$scale0.25$iterations, "iterations"))
  plot(scales[2:length(scales)], accu_nuc_total, xlab="Scale", ylab="Accuracy", main=paste("Accuracy of Nucleus by Scale",
                                                                                           masterlist$scale0.25$iterations, "iterations"))
}

# plots average cost for each iteration
plot_avgcost <- function(){
  load("~/rTest/output/1.5scale_0.25inc_10000iter.RData")
  cost_total <- NULL
  maxscale <- 1.5
  range01 <- function(x){ (x - min(x))/(max(x)-min(x)) * maxscale }
  scales <- range01(1:(length(masterlist)+1))
  for(i in masterlist){
    costs <- i$cost
    cost_total <- c(cost_total,mean(costs))
  }
  plot(scales[2:length(scales)], cost_total, xlab="Scale", ylab="Average Hungarian Cost", main=paste("Average Hungarian Cost by Scale",
                                                                                                     masterlist$scale0.25$iterations, "iterations"))
}

# Will find assignment for specific category - resid and nucleus
get_assignment_by_cat <- function(){
  x <- which(refdata$resid==1&refdata$nucleus=="C1'")
  print(refdata[x,])
  
  joined_cs <- data.frame()
  for(i in 1:nrow(predcs)){
    #print(predcs[i,2])
    ref_cs <- refdata[which(refdata$resid==as.integer(predcs[i,3])
                            &refdata$nucleus==as.character(predcs[i,2])),4]
    #print(ref_cs)
    cs_diff <- ref_cs-as.numeric(predcs[i,4])
    newrow <- data.frame(resid=as.integer(predcs[i,3]), nucleus=as.character(predcs[i,2])
                         ,cs_pred=as.numeric(predcs[i,4]), cs=ref_cs, diff=cs_diff, error=predcs[i,6],
                         std_error=cs_diff/predcs[i,6])
    if(nrow(joined_cs) == 0)
      joined_cs <- newrow
    else
      joined_cs = rbind(joined_cs, newrow)
  }
  print(joined_cs)
  return(joined_cs)
}

# plot chemical shift of predicted and measured on same graph (x,y), x = predicted, y = measured
plot_joint <- function(){
  joined_cs <- get_assignment_by_cat()
  #layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
  plot(joined_cs$cs_pred, joined_cs$cs, xlab="Predicted Chemical Shift", 
       ylab="Measured Chemical Shift", main="Predicted vs Measured Chemical Shift for Specific Residue and Nucleus")
  max_cs <- max(joined_cs$cs_pred, joined_cs$cs)
  segments(x0=0,y0=0,x1=as.integer(max_cs),y1=as.integer(max_cs))
  return(joined_cs)
}

# boxplot of each category, plotting standard errors in chemical shift 
plot_histo <- function(){
  joined_cs <- get_assignment_by_cat()
  errbox <- boxplot(joined_cs$std_error, ylab="Standard Error", main="Histogram of Standard Errors")
  print(errbox)
  print(length(errbox$out))
  print(1-length(errbox$out)/errbox$n)
}



# find most probable assignment, returns cs shift, iteration, scale
# takes in optional arguments
#      iter - numeric vector of specific iterations to look at
#      scale - numeric vector of specific scales to look at
#      ** ignores any iterations or scales that are not present
# DEFAULT: Looks at all possible iteration and scales
# Now works with specific case, 50000 iter, scale = 1, max_scale = 1
# Problem: Duplicate assignments can occur
find_probable <- function(iter=NULL,scale=NULL,max_scale=1,method='mean'){
  if(is.null(iter))
    iter <- c(10,100,1000,5000,10000)
  if(is.null(scale))
    scale <- c(0.25,0.5,0.75,1,1.25,1.5)
  
  most_probable <- list(cs=predcs, assigned=NA, iter=NA, scale=NA, assign_error=Inf)
  scale_inc <- max_scale / length(scale)
  
  for(i in iter){
    load(paste("~/rTest/output/",max_scale,"scale_",scale_inc,"inc_",i,"iter.RData",sep=""))
    for(s in scale){
      predcs_prob <- masterlist[[paste("scale",s,sep="")]]$probs
      predcs$assigned <- refdata[apply(predcs_prob, 1, which.max),"cs"]
      tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
      
      if(method=='mean'){
        abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = function(x){mean(abs(x$cs.x-x$assigned))})
        if(mean(abc$V1) < most_probable$assign_error)
          most_probable <- list(cs=predcs, assigned=abc, iter=i, scale=s, assign_error=mean(abc$V1))
      }
      else if(method=='cor'){
        abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = function(x){cor(x$cs.x,x$assigned)})
        
        if(abs(1-mean(abc$V1, na.rm=TRUE)) < most_probable$assign_error)
          most_probable <- list(cs=predcs, assigned=abc, iter=i, scale=s, assign_error=abs(1-mean(abc$V1, na.rm=TRUE)))
      }
      else if(method=='stderr'){
        abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = function(x){mean(abs((x$cs.x-x$assigned)/ x$error))})
        
        if(mean(abc$V1) < most_probable$assign_error)
          most_probable <- list(cs=predcs, assigned=abc, iter=i, scale=s, assign_error=mean(abc$V1))
      }
      
    }
  }
  return(most_probable)
}

# Given a number iterations, plot error by scale
plot_probable_byscale <- function(iter){
  scale <- c(0.25,0.5,0.75,1,1.25,1.5)
  error <- NULL
  for(s in scale){
    obj <- find_probable_hung(iter=iter,scale=s)
    error <- c(error,obj$assign_error)
  }
  plot(scale,error,main=paste("Assignment Error for",iter,"iterations of Hungarian Algorithm at each scale"))
}



'
Idea 1: same as before, but when using max, use k-largest when assignment has already been made, where k is the number of previously assigned chemical when matching the k - largest probability

i.e 3 -> 7, 4 ->7 but since 7 is already used 4 -> 9

Idea 2: Do 1 - probability matrix, apply Hungarian algorithm to find smallest cost assignment, use the returned output
'

# idea 2 Hungarian Algorithm on probabilities
# need to try 1.5
find_probable_hung <- function(iter=NULL,scale=NULL,max_scale=1.5,method='mean'){
  if(is.null(iter))
    iter <- c(10,100,1000,5000,10000)
  if(is.null(scale))
    scale <- c(0.25,0.5,0.75,1,1.25,1.5)
  
  most_probable <- list(cs=predcs, assigned=NA, iter=NA, scale=NA, assign_error=Inf)
  scale_inc <- max_scale / length(scale)
  
  for(i in iter){
    load(paste("~/rTest/output/",max_scale,"scale_",scale_inc,"inc_",i,"iter.RData",sep=""))
    for(s in scale){
      predcs_prob <- masterlist[[paste("scale",s,sep="")]]$probs
      a <- solve_LSAP(1 - predcs_prob)
      predcs$assigned <- refdata[a,"cs"]
      #predcs$assigned <- refdata[apply(predcs_prob, 1, which.max),"cs"]
      tmp <- merge(refdata, predcs, by=c("resid","nucleus"))
      
      if(method=='mean'){
        abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = function(x){mean(abs(x$cs.x-x$assigned))})
        if(mean(abc$V1) < most_probable$assign_error)
          most_probable <- list(cs=predcs, assigned=abc, iter=i, scale=s, assign_error=mean(abc$V1))
      }
      else if(method=='cor'){
        abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = function(x){cor(x$cs.x,x$assigned)})
        
        if(abs(1-mean(abc$V1, na.rm=TRUE)) < most_probable$assign_error)
          most_probable <- list(cs=predcs, assigned=abc, iter=i, scale=s, assign_error=abs(1-mean(abc$V1, na.rm=TRUE)))
      }
      else if(method=='stderr'){
        abc <- ddply(.data = tmp, .variables = c("nucleus"), .fun = function(x){mean(abs((x$cs.x-x$assigned)/ x$error))})
        
        if(mean(abc$V1) < most_probable$assign_error)
          most_probable <- list(cs=predcs, assigned=abc, iter=i, scale=s, assign_error=mean(abc$V1))
      }
      
    }
  }
  return(most_probable)
  
}

# get correlation data between assigned and actual
get_correlation <- function(method="mean", nuclei=c("H","C")){
  res <- find_probable_hung(method=method)
  # do more stuff, need to merge predcs with actual assignment and find correlation data with assigned ref and actual ref
  tmp <- merge(res$cs,refdata, by=c("resid","nucleus"),suffixes=c(".pred",".real"))
  dataCorr <- data.frame(nuc="B1",r=0,rsquare=0,RMSE=0,MAE=0,stringsAsFactors=FALSE)
  #names(dataCorr) <- c("r","r.squared","RMSE","MAE")
  for(nuc in nuclei){
    modeltmp <- lm(assigned~cs.real, data=tmp[substr(tmp$nucleus,1,1)==nuc,])
    #plot(modeltmp)
    #print(paste(n,"Nuclei"))
    #print(summary(modeltmp))
    rSquare <- summary(modeltmp)$r.squared
    if(nrow(dataCorr) == 0)
      dataCorr <- c(as.character(nuc),sqrt(rSquare),rSquare,sqrt(mean((tmp$assigned-tmp$cs.real)^2)),
                    mean(abs(tmp$assigned-tmp$cs.real)))
    else
      dataCorr = rbind(dataCorr, c(nuc,sqrt(rSquare),rSquare,sqrt(mean((tmp$assigned-tmp$cs.real)^2)),
                                 mean(abs(tmp$assigned-tmp$cs.real))))
    #print(paste("RMSE:",sqrt(mean((tmp$assigned-tmp$cs.real)^2))))
    #print(paste("MAE:",mean(abs(tmp$assigned-tmp$cs.real))))
  }
  dataCorr <- dataCorr[rownames(dataCorr) > 1,]
  return(dataCorr)
}

# plot the data in get_correlation
plot_correlation <- function(){
  #plot(modeltmp$model)
  #abline(modeltmp)
}

plot_probable_byiter <- function(scale){
  iter <- c(10,100,1000,5000,10000)
  error <- NULL
  for(i in iter){
    obj <- find_probable_hung(iter=i,scale=scale)
    error <- c(error,obj$assign_error)
  }
  plot(iter,error,main=paste("Assignment Error for",scale,"scale of Hungarian Algorithm at each iteration"))
}