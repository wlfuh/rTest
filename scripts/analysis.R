analyze <- function(modeldata, ref="actual", comp="assigned"){
  require(plyr)  
  modelerror <- ddply(.dat=modeldata, .var=c("state"), 
                      .fun = function(modeldata){
                        errors <- data.frame(state=0, mae=0, rmse=0, r=0)
                        errors$state <- modeldata$state[1]
                        # TODO MAE, RMSE
                        errors$mae <-  mean(modeldata$error*abs(modeldata[,comp]-modeldata[,ref]))
                        errors$rmse <- sqrt(mean(modeldata$error*modeldata$error*(modeldata[,comp]-modeldata[,ref])^2))
                        errors$r <- 1-cor(modeldata[,comp], modeldata[,ref], method = "pearson")
                        errors$rho <- 1-cor(modeldata[,comp], modeldata[,ref], method = "spearman")
                        errors$tau <- 1-cor(modeldata[,comp], modeldata[,ref], method = "kendall")
                        errors
                      })
  return(modelerror)
}

merge_data <- function(rna, iter, states){
	for (i in states){
	  tmp <- read.table(paste("assignments/assigned_shifts_", rna, "_", iter, "_", i, ".txt", sep =""), col.names = c("state", "resid", "resname", "nucleus", "actual", "assigned", "predicted"))
		if (i==1){
			cs <- tmp
		} else {
			cs <- rbind(cs, tmp)
		}		
	} 	
	return(cs)
}

summarize <- function(rna, iters=c(1000, 2000, 3000), states = 1:18, errortype="rmse", ref="actual", comp="assigned" ){
  cors <- rmsds <- mean_rmsd <- NULL
  for (iter in iters){
		cs <- merge_data(rna, iter, states)
		accu <- read.table("data/larmord_accuracy_resname_nucleus.txt")
		colnames(accu) <- c("resname", "nucleus", "error")
		cs <- merge(cs, accu, by = c("resname", "nucleus"))
		
		rest <- analyze(cs, ref, comp)
		struct_info <- read.table(paste("../struct_info/", rna, ".txt", sep = ""), col.names = c("ID", "decoy", "rmsd", "tmscore", "gdt", "state"))
		rest <- merge(rest, struct_info)			
		cors <- c(cors, cor(rest[,errortype], rest$rmsd, method="kendall"))
		rmsds <- c(rmsds, head(rest[order(rest[, errortype]), "rmsd"],1))
		mean_rmsd <- c(mean_rmsd, mean(head(rest[order(rest[, errortype]), "rmsd"],3)))
		#print(rest[order(rest[, errortype]),])
  }
  res <- data.frame(iters, cors, rmsds, mean_rmsd)
  return(res)
}


#summarize("1Z2J", seq(1000, 50000,1000), errortype="mae", ref = "predicted", comp = "assigned")

summarize_rnas <- function(start=1000, stop=5000, stride=1000, states=1:15, errortype="tau", rnas = unlist(strsplit("1KKA 1L1W 1OW9 1PJY 1SCL 1YSV 1Z2J 2KOC 2L3E 2LBJ 2LHP 2LI4 2LK3 2LUB 2M24 2M4W 2MFD 2MHI 2MNC 2N6S 2QH4 2Y95 5A18 5KQE", " "))){
	for (i in seq_along(rnas)){
		rna <- rnas[i]
		tmp <- summarize(rna, seq(start, stop, stride), states=states, errortype=errortype, ref = "actual", comp = "assigned")
		tmp$id <- rna
		print(tmp)
		if (i==1){
			summ <- tmp
		} else {
			summ <- rbind(summ, tmp)
		}			
	}
	return(summ)
}

