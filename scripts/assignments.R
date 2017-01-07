assign <- function(x, y){
  require(clue)
  # assigns elements in y to elements in x 
  # Input -- x (vector): e.g., chemical shifts 
  # Input -- y (vector): e.g., chemical shifts which are to be mapped to those in x (note that x can be greater than y)
  if(length(x) < length(y)){
    temp = x
    x = y
    y = temp
  }
  xmat <- matrix(x,nrow=length(y),ncol=length(x),byrow=T)
  ymat <- matrix(y,nrow=length(y),ncol=length(x),byrow=F)
  costmat <- abs(xmat-ymat)
  a <- solve_LSAP(costmat)
  
  a_mat <- matrix(0,nrow=length(y),ncol=length(x),byrow=T)
  for (i in 1:length(y)){
    j <- as.vector(a)[i]
    a_mat[i,j] <-1
  }
  return(list(a=a,a_mat=a_mat,costmat=costmat))
}
get_assignment_cost <- function(a, costmat){
  require(clue)
  # get cost of assignment
  # Input -- a (vector): linear assignment object returned by solve_LSAP
  # Input -- costmat (matrix): cost matrix
  cost = as.integer(0)
  for(i in 1:length(a))
    cost = cost + costmat[i, a[i]]
  return(cost)
}