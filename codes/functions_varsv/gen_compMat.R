gen_compMat <- function(A, n, plag){
  Jm          <- matrix(0, n*plag, n)
  Jm[1:n,1:n] <- diag(n)
  
  A   <- A[1:(n*plag),,drop=FALSE]
  Cm  <- matrix(0, n*plag, n*plag)
  if(plag==1) Cm <- t(A) else {
    for(jj in 1:(plag-1)){
      Cm[(jj*n+1):(n*(jj+1)),(n*(jj-1)+1):(jj*n)] <- diag(n)
    }
  }
  bbtemp <- A[1:(n*plag),]
  splace <- 0
  for(pp in 1:plag){
    for(nn in 1:n) {
      Cm[nn,((pp-1)*n+1):(pp*n)] <- t(bbtemp[(splace+1):(splace+n),nn])
    }
    splace <- splace+n
  }
  return(list(Cm=Cm,
              Jm=Jm))
}
