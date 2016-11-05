 library(plyr)
assign_peaks <- function(x, y){
  require(clue)
  # maps chemical shifts in y to chemical shifts in x 
  # Input -- x (vector): chemical shifts 
  # Input -- y (vector): chemical shifts which are to be mapped to those in x (not x can be greater than y)
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
    <- solve_LSAP(abs(xmat-ymat))
  cbind(x[assignments],y)
}

get_score <- function(x,y){
  result <- assign_peaks(x$predCS,y)
  cor(result[,1],result[,2],method="spearman")
  mean(abs(result[,1]-result[,2]))
}

# read in test files
csfile <- "shifts.txt"
rmsdfile <- "trajrmsd.txt"
cs <- read.table(csfile,col.names=c("processor","model","resid","resname","nucleus","predCS","expCS","ranCS","system"))
rmsd <- read.table(rmsdfile,col.names = c("model", "rmsd"))

# get scrambled chemical shifts 
actual_shifts <- subset(cs,model==1)$expCS
actual_shifts <- sample(actual_shifts)

# assign scrambled chemical shift to predicted chemical shifts and get errors
scores <- ddply(.dat=cs,.var=c("model"),.fun=get_score, y=actual_shifts)
scores$rmsd <- rmsd$rmsd
colnames(scores) <- c("model","score","rmsd")
scores <- scores[order(scores$score),]


# calculate chemical shifts from structure
