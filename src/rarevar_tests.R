
burden.fun<-function(U.sum,V.sum){

burden.pval <- pchisq(U.sum^2/V.sum, df=1, lower.tail=FALSE)
out<-data.frame(Burden.Score=U.sum, Burden.Variance=V.sum, Burden.pval=burden.pval)
class(out$Burden.Score) <- class(out$Burden.Variance) <- class(out$Burden.pval) <- "numeric"
return(out)
}

skat.fun<-function(U,V,n.variants){

Q                <- sum(U^2)
SKAT.pval        <- NA
SKAT.pval.method <- NA
if(mean(abs(V)) >= sqrt(.Machine$double.eps)) {
   pv               <- GENESIS:::.regular(Q, V, n.variants)
   SKAT.pval        <- pv$pval
   SKAT.pval.method <- pv$method
}

out<-data.frame(SKAT.pval=SKAT.pval,SKAT.pval.method=SKAT.pval.method)
class(out$SKAT.pval)        <- "numeric"
class(out$SKAT.pval.method) <- "character"
return(out)
}


skato.fun<-function(U,V,rho){

Q                 <- sum(U^2)
SKATO.pval        <- NA
SKATO.pval.method <- NA
if(mean(abs(V)) >= sqrt(.Machine$double.eps)) {
   res_skato       <- GMMAT:::.skato_pval(U = U, V = V, rho = rho, method = "davies")
   Burden.Score    <- res_skato$Burden.score
   Burden.Variance <- res_skato$Burden.var
   Burden.pval     <- res_skato$Burden.pval
   SKAT.pval       <- res_skato$SKAT.pval
   SKATO.pval      <- res_skato$p
   SKATO.minp      <- res_skato$minp
   SKATO.minp.rho  <- res_skato$minp.rho
}
out<-data.frame(Burden.Score=Burden.Score,Burden.Variance=Burden.Variance, Burden.pval=Burden.pval,SKAT.pval=SKAT.pval,SKATO.pval=SKATO.pval,SKATO.minp=SKATO.minp,SKATO.minp.rho=SKATO.minp.rho)
return(out)
}


smmat.fun<-function(U,V,U.sum,V.sum,GG1){
# Compute burden-adjusted SKAT statistic
  U <- U - GG1*U.sum/V.sum
  Q <- sum(U^2)
  V <- V - tcrossprod(GG1)/V.sum
# SKAT
  theta.pval        <- NA
  theta.pval.method <- NA
  err               <- NA
  if(mean(abs(V)) >= sqrt(.Machine$double.eps)) {
     pv                <- GENESIS:::.regular(Q, V, n.variants)
     theta.pval        <- pv$pval
     theta.pval.method <- pv$method
     err               <- pv$err
  }
# Fisher's method to combine p-values
  SMMAT.pval <- tryCatch(pchisq(-2*log(burden.pval)-2*log(theta.pval), df=4, lower.tail = FALSE),
                         error = function(e) { NA })
  if(is.na(SMMAT.pval)) {
     err        <- 1
     SMMAT.pval <- NA
     SMMAT.pval <- burden.pval
  }
  out<-data.frame(theta.pval=theta.pval, theta.pval.method=theta.pval,err=err, SMMAT.pval=SMMAT.pval)
  class(out$theta.pval) <- class(out$err) <- class(out$SMMAT.pval) <- "numeric"
  class(out$theta.pval.method) <- "character"

return(out)
}
