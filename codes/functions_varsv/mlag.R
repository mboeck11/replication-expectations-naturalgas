mlag <- function(Y,plag){
  Y <- as.matrix(Y)
  Traw <- nrow(Y)
  n <- ncol(Y)
  Ylag <- matrix(0,Traw,plag*n)
  for (ii in 1:plag){
    Ylag[(plag+1):Traw,(n*(ii-1)+1):(n*ii)] <- Y[(plag+1-ii):(Traw-ii),(1:n)]
  }
  colnames(Ylag) <- paste0(colnames(Y),".lag",rep(seq(plag),each=n))
  return(Ylag)
}
