pathname <- "/home/afrankz/testbed/testWilliam/calcs/"
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

merge_data <- function(rna, predictor, iter, states){
	for (i in states){
	  tmp <- read.table(paste(pathname,"assignments/assigned_shifts_", predictor, "_", rna, "_", iter, "_", i, ".txt", sep =""), col.names = c("state", "resid", "resname", "nucleus", "actual", "assigned", "predicted"))
		if (i==1){
			cs <- tmp
		} else {
			cs <- rbind(cs, tmp)
		}		
	} 	
	return(cs)
}

summarize <- function(rna, predictor, iters=c(1000, 2000, 3000), states = 1:18, errortype="rmse", ref="actual", comp="assigned" ){
  cors <- rmsds <- mean_rmsd <- NULL
  accu_file <- paste(pathname,"/data/", predictor,"_accuracy_resname_nucleus.txt", sep="")
  for (iter in iters){
		cs <- merge_data(rna, predictor, iter, states)
		accu <- read.table(accu_file)
		colnames(accu) <- c("resname", "nucleus", "error")
		cs <- merge(cs, accu, by = c("resname", "nucleus"))
		
		rest <- analyze(cs, ref, comp)
		struct_info <- read.table(paste(pathname,"../struct_info/", rna, ".txt", sep = ""), col.names = c("ID", "decoy", "rmsd", "tmscore", "gdt", "state"))
		rest <- merge(rest, struct_info)			
		cors <- c(cors, cor(rest[,errortype], rest$rmsd, method="kendall"))
		rmsds <- c(rmsds, head(rest[order(rest[, errortype]), "rmsd"],1))
		mean_rmsd <- c(mean_rmsd, mean(head(rest[order(rest[, errortype]), "rmsd"],3)))
		#print(rest[order(rest[, errortype]),])
  }
  res <- data.frame(iters, cors, rmsds, mean_rmsd)
  return(res)
}


summarize_rnas <- function(ref="actual",predictor="ramsey", start=1000, stop=10000, stride=1000, states=1:15, errortype="tau", rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2N6W 2NCI 1UUU 2JWV 2M12 2M8K 2N2O 2N6X 2QH2 5A17 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 5A18 2M22 2MFD 2N4L 5KQE 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 1R7W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2LU0 2Y95 1SCL 1Z2J 2L3E 2M24 2MHI 2N6S", " "))){
	for (i in seq_along(rnas)){
		rna <- rnas[i]
		tmp <- summarize(rna, predictor, seq(start, stop, stride), states=states, errortype=errortype, ref = ref, comp = "assigned")
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
ref="predicted"
methods=c("rho","tau","r","mae","rmse")
for(method in methods){
  results <- summarize_rnas(ref=ref,rnas = unlist(strsplit("1ZC5 2LUN 2M4W 2MNC 2N6T 2M5U 2MXL 2N6W 2NCI 1UUU 2JWV 2M12 2M8K 2N2O 2N6X 2QH2 5A17 1XHP 2K66 2LQZ 2M21 2MEQ 2N2P 2N7X 2QH4 5A18 2M22 2MFD 2N4L 5KQE 1KKA 1PJY 2LI4 2LUB 4A4S 1L1W 1R7W 2LBJ 2LK3 4A4T 1LC6 2FDT 2LBL 2LP9 2LV0 4A4U 1LDZ 2LDL 2LPA 1NC0 2LDT 1OW9 1YSV 2KOC 2LHP 2Y95 1SCL 1Z2J 2L3E 2M24 2MHI 2N6S"," "))
  ,errortype=method)
  ramsey_training_set <- unlist((strsplit("1LDZ 1LC6 1NC0 2KOC 1PJY 1OW9 1R7W 1R7Z 1YSV 1Z2J 2GM0 2JTP 2JXQ 2JXS 2K3Z 2K41 2L3E 2LDL 1UUU", " ")))
  training_results <- results[(results$id %in% ramsey_training_set),]
  testing_results <- results[!(results$id %in% ramsey_training_set),]
  final <- subset(testing_results,iters==10000)
  final$status <- final$rmsds < 3.0
  finalall <- results
  finalall$status <- finalall$rmsds < 3.0
  mean(final$status)
  writedir=paste("data/",method,"/",sep="")
  write.table(final,file=paste(writedir,ref,"_results_final_10000.txt",sep=""),quote=FALSE)
  write.table(finalall,file=paste(writedir,ref,"_results_final_all.txt",sep=""),quote=FALSE)
  write.table(training_results,file=paste(writedir,"training_results_",ref,".txt",sep=""),quote=FALSE)
  write.table(testing_results,file=paste(writedir,"testing_results_",ref,".txt",sep=""),quote=FALSE)
}
# how many if (nrow(unique(read.table("data/observed_shifts_1SCL.txt"))) < nrow(unique(read.table("data/observed_shifts_1SCL.txt")))

