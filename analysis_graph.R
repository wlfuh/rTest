# Things to look into

# Assigned vs Actual  and Assigned vs Predicted

# Things to compare
# Impact of Iterations
# Iteration vs correlation of rmsd
# Iteration vs mean rmsd value

# RMSD Correlation (10000 iterations, tau method)
# Mean RMSD Value
# Correlation error value
# Histogram with each ID
# Summary of top ID's with highest cors or mean_rmsd

# By Iterations (tau method)
# Over all id's values averaged
# Two graphs of correlation and mean rmsd value

# By error method (10000 iteration)
# Over all id's values averaged
# Two graphs of correlation and mean rmsd value

require(ggplot2)
require(plyr)

plot_iter <- function(result, field="cors", doplot=TRUE){
  result <- result[order(result[,"iters"]),]
  data <- NULL
  iterations <- unique(result$iters)
  for(iter in iterations){
    data <- c(data, mean(result[result$iters==iter, field]))
  }
  if(field=="cors" && doplot)
    plot(iterations, data, xlab="Iterations", main="RMSD Correlation over Iterations", ylab = "RMSD Correlation")
  else if(field=="mean_rmsd" && doplot)
    plot(iterations, data, xlab="Iterations", main="Mean RMSD over Iterations", ylab = "Mean RMSD (10e-10 m)")
  return(data.frame(iterations, data))
}



allRNA <- function(type="actual"){
  errortypes <- c("mae","r","rho","rmse","tau")
  iterations <- c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
  cor_combined <- NULL
  mrmsd_combined <- NULL
  variables <- NULL
  # do graph by iterations for each error
  for(errortype in errortypes){
    res <- read.table(file=paste("data/",errortype,"/",type,"_results_final_all.txt",sep=""))
    png(paste("graphs/assigned_",type,"/",errortype,"/","iteration_rmsdcorr_graph.png",sep=""))
    data <- plot_iter(res, field="cors")
    cor_combined <- c(cor_combined, data$data)
    dev.off()
    write.table(data,file=paste("graphs/assigned_",type,"/",errortype,"/","iteration_rmsdcorr_data.txt",sep=""),row.names=FALSE,quote=FALSE)
    png(paste("graphs/assigned_",type,"/",errortype,"/","iteration_rmsdmean_graph.png",sep=""))
    data <- plot_iter(res, field="mean_rmsd")
    mrmsd_combined <- c(mrmsd_combined, data$data)
    dev.off()
    write.table(data,file=paste("graphs/assigned_",type,"/",errortype,"/","iteration_rmsdmean_data.txt",sep=""),row.names=FALSE,quote=FALSE)
    variables <- c(variables, rep(errortype, length(iterations)))  
  }
  
  cor_graph <- rename(data.frame(rep(iterations, length(errortypes)), cor_combined, variables)
                      , c("rep.iterations..length.errortypes.."="iterations"))
  ggplot(data = cor_graph, aes(x=iterations, y=cor_combined)) + geom_line(aes(colour=variables)) + xlab("Iterations") + ylab("RMSD Correlation") + ggtitle("RMSD Correlation over Iterations")
  ggsave(paste("graphs/assigned_",type,"/combined/iteration_rmsdcorr_graph.png",sep=""),
         width = 6.7, height = 4.16)
  write.table(cor_graph
              ,file=paste("graphs/assigned_",type,"/combined/iteration_rmsdcorr_data.txt",sep=""),row.names=FALSE,quote=FALSE)
  
  mrmsd_graph <- rename(data.frame(rep(iterations, length(errortypes)), mrmsd_combined, variables)
                      , c("rep.iterations..length.errortypes.."="iterations"))
  ggplot(data = mrmsd_graph, aes(x=iterations, y=mrmsd_combined)) + geom_line(aes(colour=variables)) + xlab("Iterations") + ylab("Mean RMSD") + ggtitle("Mean RMSD over Iterations")
  ggsave(paste("graphs/assigned_",type,"/combined/iteration_rmsdmean_graph.png",sep=""),
         width = 6.7, height = 4.16)
  write.table(mrmsd_graph
              ,file=paste("graphs/assigned_",type,"/combined/iteration_rmsdmean_data.txt",sep=""),row.names=FALSE,quote=FALSE)
}

randomRNA <- function(type="actual"){
  RNA_seq <- c('2N7X', '2N6T', '1R7W', '4A4U', '1LDZ', '2N6S', '2QH4', '1NC0', '1UUU', '2LQZ')
  errortypes <- c("mae","r","rho","rmse","tau")
  iterations <- c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
  for(rna in RNA_seq){
    cor_combined <- NULL
    mrmsd_combined <- NULL
    variables <- NULL
    for(errortype in errortypes){
      res <- read.table(file=paste("data/",errortype,"/",type,"_results_final_all.txt",sep=""))
      data <- plot_iter(res[res$id==rna,], field="cors", doplot=FALSE)
      cor_combined <- c(cor_combined, data$data)
      data <- plot_iter(res[res$id==rna,], field="mean_rmsd", doplot=FALSE)
      mrmsd_combined <- c(mrmsd_combined, data$data)
      variables <- c(variables, rep(errortype, length(iterations)))
    }
    cor_graph <- rename(data.frame(rep(iterations, length(errortypes)), cor_combined, variables)
                        , c("rep.iterations..length.errortypes.."="iterations"))
    ggplot(data = cor_graph, aes(x=iterations, y=cor_combined)) + geom_line(aes(colour=variables)) + xlab("Iterations") + ylab("RMSD Correlation") + ggtitle(paste("RMSD Correlation over Iterations for", rna))
    ggsave(paste("graphs/assigned_",type,"/combined/",rna,"_iteration_rmsdcorr_graph.png",sep=""),
           width = 6.7, height = 4.16)
    write.table(cor_graph
                ,file=paste("graphs/assigned_",type,"/combined/",rna,"_iteration_rmsdcorr_data.txt",sep=""),row.names=FALSE,quote=FALSE)
    
    mrmsd_graph <- rename(data.frame(rep(iterations, length(errortypes)), mrmsd_combined, variables)
                          , c("rep.iterations..length.errortypes.."="iterations"))
    ggplot(data = mrmsd_graph, aes(x=iterations, y=mrmsd_combined)) + geom_line(aes(colour=variables)) + xlab("Iterations") + ylab("Mean RMSD") + ggtitle(paste("Mean RMSD over Iterations for", rna))
    ggsave(paste("graphs/assigned_",type,"/combined/",rna,"_iteration_rmsdmean_graph.png",sep=""),
           width = 6.7, height = 4.16)
    write.table(mrmsd_graph
                ,file=paste("graphs/assigned_",type,"/combined/",rna,"_iteration_rmsdmean_data.txt",sep=""),row.names=FALSE,quote=FALSE)
  }
  
}

measure_acc <- function(lim=2, type="actual"){
  # create dataframe with % accuracy for each error method
  # also write id of "acceptable" rmsd
  errortypes <- c("mae","r","rho","rmse","tau")
  accuracy <- NULL
  for(errortype in errortypes){
    res <- read.table(file=paste("data/",errortype,"/",type,"_results_final_10000.txt",sep=""))
    sub_res <- res[res$rmsds <= lim,]
    accuracy <- c(accuracy, nrow(sub_res)/nrow(res))
    write.table(sub_res,file=paste("analysis/assigned_",type,"/",errortype,"/","acceptable_rna_lim_",lim,".txt",sep=""),row.names=FALSE,quote=FALSE)
  }
  acc <- data.frame(errortypes, accuracy)
  write.table(acc,file=paste("analysis/assigned_",type,"/accuracy_lim_",lim,".txt",sep=""),row.names=FALSE,quote=FALSE)
  return(acc)
}

plot_acc <- function(type="actual"){
  steps <- seq(0, 5, 0.5)
  accuracy <- NULL
  errortype <- NULL
  limit <- NULL
  for(s in steps){
    dat <- measure_acc(lim=s, type)
    accuracy <- c(accuracy, dat$accuracy)
    errortype <- c(errortype, as.vector(dat$errortypes))
    limit <- c(limit, rep(s, nrow(dat)))
  }
  acc_graph <- data.frame(limit, errortype, accuracy)
  ggplot(data = acc_graph, aes(x=limit, y=accuracy)) + geom_line(aes(colour=errortype)) + xlab("RMSD Limit") + ylab("Accuracy (0-1)") + ggtitle(paste("Accuracy of assigned vs ", type, "chemical shift based on RMSDS Limit"))
  ggsave(paste("analysis/assigned_",type,"/graphs/",type,"_accuracy_rmsds_graph.png",sep=""),
         width = 6.7, height = 4.16)
  write.table(acc_graph
              ,file=paste("analysis/assigned_",type,"/graphs/",type,"_accuracy_rmsds_data.txt",sep=""),row.names=FALSE,quote=FALSE)
}

run_tests <- function(){
  types <- c("predicted", "actual")
  lims <- c(2,2.5,3)
  for(type in types){
    #allRNA(type=type)
    #randomRNA(type=type)
    #for(l in lims)
    #  measure_acc(lim=l,type=type)
    plot_acc(type)
  }
}
