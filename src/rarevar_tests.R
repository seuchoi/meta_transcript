
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
   pv               <- regular(Q, V, n.variants)
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
   res_skato       <- skato_pval(U = U, V = V, rho = rho, method = "davies")
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


smmat.fun<-function(U,V,U.sum,V.sum,GG1,n.variants){
# Compute burden-adjusted SKAT statistic
  U <- U - GG1*U.sum/V.sum
  Q <- sum(U^2)
  V <- V - tcrossprod(GG1)/V.sum
  burden.pval <- pchisq(U.sum^2/V.sum, df=1, lower.tail=FALSE)

  # SKAT
  theta.pval        <- NA
  theta.pval.method <- NA
  err               <- NA
  if(mean(abs(V)) >= sqrt(.Machine$double.eps)) {
     pv                <- regular(Q, V, n.variants)
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





regular<-function (Q, V, ncolG) { ### From GENESIS
    if (ncolG == 1) {
        pv <- list(pval = pchisq(as.numeric(Q/V), df = 1, lower.tail = FALSE), 
            method = "integration")
    }
    else {
        lambda <- eigen(V, only.values = TRUE, symmetric = TRUE)$values
        pv <- .pchisqsum(x = Q, df = rep(1, length(lambda)), 
            a = lambda)
    }
    pv$err <- ifelse(is.na(pv$pval), 1, 0)
    return(pv)
}

	   
Q_pval <- function(Q, lambda, method = "davies") {
  if(method == "davies") {
    tmp <- try(suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6)))
    if(inherits(tmp, "try-error") || tmp$ifault > 0 || tmp$Qq <= 1e-5 || tmp$Qq >= 1) method <- "kuonen"
    else return(tmp$Qq)
  }
  if(method == "kuonen") {
    pval <- try(.pKuonen(x = Q, lambda = lambda))
    if(inherits(pval, "try-error") || is.na(pval)) method <- "liu"
    else return(pval)
  }
  if(method == "liu") {
    pval <- try(CompQuadForm::liu(q = Q, lambda = lambda))
    if(inherits(pval, "try-error")) cat("Warning: method \"liu\" failed...\nQ:", Q, "\nlambda:", lambda, "\n")
    else return(pval)
  }
  return(NA)
}


skato_pval <- function(U, V, rho, method = "davies") { ## from GMMAT
    n.r <- length(rho)
    n.p <- length(U)
    lambdas <- vector("list", n.r)
    pval <- qval <- rep(NA, n.r)
    Q <- (1-rho)*sum(U^2)+rho*sum(U)^2
    Burden.score <- Burden.var <- Burden.pval <- SKAT.pval <- NA
    for(i in 1:n.r) {
	if(rho[i]==1) {
	    Burden.score <- sum(U)
	    Burden.var <- sum(V)
	    Burden.pval <- pchisq(Burden.score^2/Burden.var, df=1, lower.tail=FALSE)
	    lambdas[[i]] <- Burden.var
	    pval[i] <- Burden.pval
	    next
	}
	if(rho[i]!=0) {
	    R.M <- matrix(rho[i], n.p, n.p)
	    diag(R.M) <- 1
	    R.M.chol <- t(chol(R.M, pivot = TRUE))
	    V.temp <- crossprod(R.M.chol, crossprod(V, R.M.chol))
	} else V.temp <- V
	lambda <- eigen(V.temp, only.values = TRUE, symmetric=TRUE)$values
    	lambdas[[i]] <- lambda[lambda > 0]
	pval[i] <- Q_pval(Q[i], lambdas[[i]], method = method)
	if(rho[i]==0) SKAT.pval <- pval[i]
    }
    minp <- min(pval)
  if(any(is.na(pval))){
       return(list(p=NA, minp=NA, minp.rho=NA, Burden.score=Burden.score, Burden.var=Burden.var, Burden.pval=Burden.pval, SKAT.pval=SKAT.pval))
    }
    for(i in 1:n.r) {
	df <- sum(lambdas[[i]]^2)^2/sum(lambdas[[i]]^4)
	qval[i] <- (qchisq(minp, df, lower.tail = FALSE)-df)/sqrt(2*df)*sqrt(2*sum(lambdas[[i]]^2))+sum(lambdas[[i]])
    }
    ZMZ <- tcrossprod(rowSums(V))/sum(V)
    V.temp <- V - ZMZ
    lambda <- eigen(V.temp, only.values = TRUE, symmetric = TRUE)$values
    lambda <- lambda[lambda > 0]
    muq <- sum(lambda)
    varq <- sum(lambda^2) * 2 + sum(ZMZ * V.temp) * 4
    df <- sum(lambda^2)^2/sum(lambda^4)
    tau <- rho * sum(V) + sum(V %*% V)/sum(V) * (1 - rho)
    re <- tryCatch({
        integrate(function(x){
    	    t1 <- tau %x% t(x)
    	    re<-pchisq((apply((qval - t1)/(1-rho),2,min) - muq)/sqrt(varq)*sqrt(2*df) + df, df=df) * dchisq(x,df=1)
    	    return(re)
	}, lower = 0, upper = 40, subdivisions = 2000, abs.tol = 10^-25)
    }, error=function(e) NA)
    return(list(p = min(1-re[[1]], minp*n.r), minp = minp, minp.rho = rho[which.min(pval)],
    Burden.score=Burden.score, Burden.var=Burden.var, Burden.pval=Burden.pval,
    SKAT.pval=SKAT.pval))
}
