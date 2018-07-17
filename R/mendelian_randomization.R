# scripts for Mendelian Randomization of stickleback armor
# original in edaII.R, started Aug 9, 2013
# April 6, 2018

mv_selection_differential <- function(Xa, Xb, codes=NULL, freq_a=NULL, freq_b=NULL, adjust=FALSE, bootstrap=FALSE, niter=4999){
  # this is a function specific to Rennison et al, where one column is genotype and the expected frequency of each level of genotype is specified instead of computed from X. And, if adjust=TRUE there is a column "length" of standard length.
  # Xa contains the data before selection, one column must be genotype with levels coded as integer
  # Xb contains the data after selection, one column must be genotype with levels coded as integer
  # codes is the coding of genotype, in the order given in freq
  # freq is a vector containing the expected frequencies of the unique values of genotype, in the order given in codes. The procedure will sample this so each iteration, the realized frequency of the samples will differ from the expectation.
  # adjust, if TRUE then condition on column "length" prior to computing selection differential (as opposed to including in linear model)
  
  n_codes <- lengths(codes)
  inc_a <- list()
  inc_b <- list()
  for(j in codes){
    inc_a[[as.character(j)]] <- which(Xa[,genotype]==j)
    inc_b[[as.character(j)]] <- which(Xb[,genotype]==j)
  }
  n_a <- nrow(Xa)
  n_b <- nrow(Xb)
  
  if(bootstrap==FALSE){niter <- 0}
  B <- B_prime <- matrix(NA, nrow=(niter+1), ncol=2) # B is the matrix of b vectors for each iteration of the bootstrap
  colnames(B) <- c('adj.plates', 'genotype')
  colnames(B_prime) <- paste(colnames(B),"_prime", sep="")
  
  for(iter in 1:(niter+1)){
    # resample inc_A and inc_B
    freq_a_i <- table(sample(codes, size=n_a, prob=freq_a, replace=TRUE))
    freq_b_i <- table(sample(codes, size=n_b, prob=freq_b, replace=TRUE))
    inc_a_i <- inc_b_i <- NULL
    for(j in codes){
      inc_a_i <- c(inc_a_i, sample(inc_a[[as.character(j)]], size=freq_a_i[as.character(j)], replace=TRUE))
      inc_b_i <- c(inc_b_i, sample(inc_b[[as.character(j)]], size=freq_b_i[as.character(j)], replace=TRUE))
    }
    
    # if size adjust
    if(adjust==TRUE){
      dt <- rbind(Xa[inc_a_i,], Xb[inc_b_i,])
      fit <- lm(raw.logplates ~ length*genotype, data=dt[length < 34.639958,])
      dt[, adj.logplates:=ifelse(length < 34,
                          predict(fit, data.frame(length=34.0, genotype=genotype)) + 
                          raw.logplates - predict(fit, data.frame(length=length, genotype=genotype)),
                          raw.logplates)]
      # limit adjusted plates to max of log(32)
      dt[, adj.logplates:=ifelse(adj.logplates > log(32), log(32), adj.logplates)]
      dt[, adj.plates:=exp(adj.logplates)]
      Xai <- dt[1:nrow(Xa), .SD, .SDcols=c('adj.plates', 'genotype')]
      Xbi <- dt[(nrow(Xa)+1):nrow(dt), .SD, .SDcols=c('adj.plates', 'genotype')]
    }else{
      Xai <- Xa[inc_a_i,]
      Xbi <- Xb[inc_b_i,]
    }
    
    S <- cov(Xai)
    xbar <- apply(Xbi, 2, mean) - apply(Xai, 2, mean)
    B[iter,] <- solve(S)%*%xbar
    sd <- sqrt(diag(S))
    B_prime[iter,] <- B[iter,]*sd
  }
  return(cbind(B, B_prime))
}

mr_selection_differential <- function(Xa, Xb, codes=NULL, freq_a=NULL, freq_b=NULL, bootstrap=FALSE, niter=4999){
  # this is a function specific to Rennison et al, where one column is genotype and the expected frequency of each level of genotype is specified instead of computed from X
  # the output is the mendelian randomization
  # Xa contains the data before selection, one column must be genotype (=Z) with levels coded as integer
  # Xb contains the data after selection, one column must be genotype with levels coded as integer
  # codes is the coding of genotype, in the order given in freq
  # freq is a vector containing the expected frequencies of the unique values of genotype, in the order given in codes. The procedure will sample this so each iteration, the realized frequency of the samples will differ from the expectation.
  xcol <- setdiff(colnames(Xa), 'genotype')
  xa <- Xa[, get(xcol)]
  xb <- Xb[, get(xcol)]
  za <- Xa[, genotype]
  zb <- Xb[, genotype]
  
  n_codes <- lengths(codes)
  inc_a <- list()
  inc_b <- list()
  for(j in codes){
    inc_a[[as.character(j)]] <- which(Xa[,genotype]==j)
    inc_b[[as.character(j)]] <- which(Xb[,genotype]==j)
  }
  n_a <- nrow(Xa)
  n_b <- nrow(Xb)
  
  if(bootstrap==FALSE){niter <- 0}
  B <- matrix(NA, nrow=(niter+1), ncol=4)
  B_prime <- matrix(NA, nrow=(niter+1), ncol=4)
  colnames(B) <- c('beta.xy', 'beta.zy', 'beta.zx', 'beta.xy_MR')
  colnames(B_prime) <- paste(colnames(B),"_prime", sep="")

  for(iter in 1:(niter+1)){
    # resample inc_A and inc_B
    freq_a_i <- table(sample(codes, size=n_a, prob=freq_a, replace=TRUE))
    freq_b_i <- table(sample(codes, size=n_b, prob=freq_b, replace=TRUE))
    inc_a_i <- inc_b_i <- NULL
    for(j in codes){
      inc_a_i <- c(inc_a_i, sample(inc_a[[as.character(j)]], size=freq_a_i[as.character(j)], replace=TRUE))
      inc_b_i <- c(inc_b_i, sample(inc_b[[as.character(j)]], size=freq_b_i[as.character(j)], replace=TRUE))
    }
    
    beta.xy <- mean(xb[inc_b_i],na.rm=TRUE) - mean(xa[inc_a_i],na.rm=TRUE)
    beta.zy <- mean(zb[inc_b_i],na.rm=TRUE) - mean(za[inc_a_i],na.rm=TRUE)
    beta.zx <- coefficients(lm(xa[inc_a_i]~za[inc_a_i]))['za[inc_a_i]']
    # beta.zx <- mean(xa[inc_a_i][za[inc_a_i]==1], na.rm=TRUE) - mean(xa[inc_a_i][za[inc_a_i]==0], na.rm=TRUE)
    mr <- beta.zy/beta.zx
    B[iter, ] <- c(beta.xy, beta.zy, beta.zx, mr)
    
    beta.xy_prime <- beta.xy/sd(xa[inc_a_i],na.rm=TRUE)
    beta.zy_prime <- beta.zy/sd(za[inc_a_i],na.rm=TRUE)
    beta.zx_prime <- beta.zx/sd(xa[inc_a_i],na.rm=TRUE)
    mr_prime <- mr/sd(xa[inc_a_i],na.rm=TRUE)
    B_prime[iter, ] <- c(beta.xy_prime, beta.zy_prime, beta.zx_prime, mr_prime)
  }

  return(cbind(B, B_prime))
}

mr_selection_differential_gen <- function(Xa, Xb, codes=NULL, freq_a=NULL, freq_b=NULL, bootstrap=FALSE, niter=4999){
  # generalization of the function built for Rennison et al
  # the output is the mendelian randomization
  # Xa contains the data before selection, one column must be genotype (=Z) with levels coded as integer
  # Xb contains the data after selection, one column must be genotype with levels coded as integer
  # codes is the coding of genotype, in the order given in freq
  # freq is a vector containing the expected frequencies of the unique values of genotype, in the order given in codes. The procedure will sample this so each iteration, the realized frequency of the samples will differ from the expectation.
  
  xcol <- setdiff(colnames(Xa), 'genotype')
  xa <- Xa[, get(xcol)]
  xb <- Xb[, get(xcol)]
  za <- Xa[, genotype]
  zb <- Xb[, genotype]
  
  n_codes <- lengths(codes)
  inc_a <- list()
  inc_b <- list()
  for(j in codes){
    inc_a[[as.character(j)]] <- which(Xa[,genotype]==j)
    inc_b[[as.character(j)]] <- which(Xb[,genotype]==j)
  }
  n_a <- nrow(Xa)
  n_b <- nrow(Xb)
  
  if(bootstrap==FALSE){niter <- 0}
  B <- matrix(NA, nrow=(niter+1), ncol=4)
  B_prime <- matrix(NA, nrow=(niter+1), ncol=4)
  colnames(B) <- c('beta.xy', 'beta.zy', 'beta.zx', 'beta.xy_MR')
  colnames(B_prime) <- paste(colnames(B),"_prime", sep="")
  
  for(iter in 1:(niter+1)){
    # resample inc_A and inc_B
    freq_a_i <- table(sample(codes, size=n_a, prob=freq_a, replace=TRUE))
    freq_b_i <- table(sample(codes, size=n_b, prob=freq_b, replace=TRUE))
    inc_a_i <- inc_b_i <- NULL
    for(j in codes){
      inc_a_i <- c(inc_a_i, sample(inc_a[[as.character(j)]], size=freq_a_i[as.character(j)], replace=TRUE))
      inc_b_i <- c(inc_b_i, sample(inc_b[[as.character(j)]], size=freq_b_i[as.character(j)], replace=TRUE))
    }
    
    beta.xy <- mean(xb[inc_b_i],na.rm=TRUE) - mean(xa[inc_a_i],na.rm=TRUE)
    beta.zy <- mean(zb[inc_b_i],na.rm=TRUE) - mean(za[inc_a_i],na.rm=TRUE)
    beta.zx <- coefficients(lm(xa[inc_a_i]~za[inc_a_i]))['za[inc_a_i]']
    # beta.zx <- mean(xa[inc_a_i][za[inc_a_i]==1], na.rm=TRUE) - mean(xa[inc_a_i][za[inc_a_i]==0], na.rm=TRUE)
    mr <- beta.zy/beta.zx
    B[iter, ] <- c(beta.xy, beta.zy, beta.zx, mr)
    
    beta.xy_prime <- beta.xy/sd(xa[inc_a_i],na.rm=TRUE)
    beta.zy_prime <- beta.zy/sd(za[inc_a_i],na.rm=TRUE)
    beta.zx_prime <- beta.zx/sd(xa[inc_a_i],na.rm=TRUE)
    mr_prime <- mr/sd(xa[inc_a_i],na.rm=TRUE)
    B_prime[iter, ] <- c(beta.xy_prime, beta.zy_prime, beta.zx_prime, mr_prime)
  }
  
  return(cbind(B, B_prime))
}


selection_differential <- function(x,y,std=FALSE){
  # x is a continuous morphological trait
  # y is NOT survival but experimental treatment, with 1=treated, 0=control
  # if x is standardized then this returns the standardized selection differential assume sd(survival) = 0.5
  # or if std=TRUE this will return standardized coefficient based on sd in untreaded group
  
  if(std==TRUE){x <- x/sd(x[y==0],na.rm=TRUE)}
  N <- length(x)
  beta.xy <- mean(x[y==1],na.rm=TRUE) - mean(x[y==0],na.rm=TRUE)
  return(beta.xy)	
}

mendelian_randomization <- function(x, y, z, std=FALSE, bootstrap=FALSE, niter=4999, A=1, B=2){
  # the original code from 2013
  # computes MR of a selection differential, B - A, Not regression!
  # B should be the selected group, A should be the control
  # x is a continuous morphological trait
  # y is NOT survival but experimental treatment
  # z is genotype score
  # A and B are group identification. The differential is B - A
  # if bootstrap then return bootstrap distribution of regression
  groups <- levels(y)
  dt <- data.table(x=x, y=y, z=z)
  A_cases <- which(y==groups[A]) # control group cases
  B_cases <- which(y==groups[B]) # selected group cases
  inc_A <- A_cases
  inc_B <- B_cases
  if(bootstrap==FALSE){niter <- 0}
  beta.xy <- numeric(niter+1)
  beta.zy <- numeric(niter+1)
  beta.zx <- numeric(niter+1)
  mr <- numeric(niter+1)
  for(iter in 1:(1+niter)){
    beta.xy[iter] <- mean(x[inc_B],na.rm=TRUE) - mean(x[inc_A],na.rm=TRUE)
    beta.zy[iter] <- mean(z[inc_B],na.rm=TRUE) - mean(z[inc_A],na.rm=TRUE)
    beta.zx[iter] <- summary(lm(x[inc_A]~z[inc_A]))$coefficients[2]
    # beta.zx <- mean(x[y==groups[A] & z==1], na.rm=TRUE) - mean(x[y==groups[A] & z==0], na.rm=TRUE)
    mr[iter] <- beta.zy[iter]/beta.zx[iter]
    if(std==TRUE){
      beta.xy[iter] <- beta.xy[iter]/sd(x[inc_A],na.rm=TRUE)
      beta.zy[iter] <- beta.zy[iter]/sd(z[inc_A],na.rm=TRUE)
      beta.zx[iter] <- beta.zx[iter]/sd(x[inc_A],na.rm=TRUE)
      mr[iter] <- mr[iter]/sd(x[inc_A],na.rm=TRUE)
    }
    inc_A <- sample(A_cases, replace=TRUE)
    inc_B <- sample(B_cases, replace=TRUE)
  }
  out <- data.frame(beta.xy=beta.xy,beta.zy=beta.zy,beta.zx=beta.zx,mr=mr)
  
  return(out)
}

size_corrected <- function(x, sl, origin=NULL){
  # size correction via residual of regression of x on standard length sl.
  # if origin is null then returns residuals, else returns predicted value at x=origin
  # to correct within groups, simply use by=.(grouping) in the data.table
  # e.g.
  # marchinko[, ant.dorspine.s:=size_corrected(ant.dorspine, standlength), by=family]
  
  fit <- lm(x ~ sl, na.action = na.exclude)
  x.adj <- residuals(fit)
  if(!is.null(origin)){
    add <- predict(fit, newdata=data.frame(sl=origin))
    x.adj <- x.adj + add
  }
  return(x.adj)
}
