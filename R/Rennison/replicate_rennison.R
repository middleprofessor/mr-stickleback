# replicate Rennison et al. 2015
# Jeffrey A. Walker
# April 8, 2018
# to send to Rennison et al. to get code - copied from armor_mr.Rmd

# libraries
library(data.table)

# my version of the data has two changes to the header row of "Rennison et al archive.csv"
# 1) there are two columns with label "genotype" in the original - I've changed the second column (with -1, 0, 1 scores) to "eda_add"
# 2) I changed "dominance" to "eda_dom"

  # open data
  file_path <- "data/Rennison-am_nat-2016.txt"
  rennison <- data.table(fread(file_path))
  rennison[, month:=factor(month, c('06-Sep', '06-Oct', '06-Nov'))]
  rennison[, genotype:=factor(genotype)]

# size adjustment doesn't replicate
# adjust using ANCOVA model including interaction
  fit <- lm(raw.logplates ~ length*genotype, data=rennison[length<34,])
  
  # if length < 34.0 mm then value is predicted value at length=34 (conditional on genotype) plus residual value (conditional on length*genotype)
  rennison[, logplates.s2:=ifelse(length < 34,
                                  predict(fit, data.frame(length=34.0, genotype=genotype)) + 
                                  raw.logplates - predict(fit, data.frame(length=length, genotype=genotype)),
                                  raw.logplates)]
  # limit adjusted plates to max of log(32)
  rennison[, logplates.s2:=ifelse(logplates.s2 > log(32), log(32), logplates.s2)]
  
# conditional effects don't replicate even when using original adjusted plates data
# I get .29
  Xa <- rennison[month==levels(month)[1] & eda_add!=-1, .SD, .SDcols=c('adj.plates', 'eda_dom')]
  Xb <- rennison[month==levels(month)[2] & eda_add!=-1, .SD, .SDcols=c('adj.plates', 'eda_dom')]
  res1 <- mv_selection_differential(Xa, Xb, std=TRUE, bootstrap=FALSE)
  # apply(res1,2, quantile, c(0.025, 0.5, 0.975))
  
  Xa <- rennison[month==levels(month)[2] & eda_add!=-1, .SD, .SDcols=c('adj.plates', 'eda_dom')]
  Xb <- rennison[month==levels(month)[3] & eda_add!=-1, .SD, .SDcols=c('adj.plates', 'eda_dom')]
  res2 <- mv_selection_differential(Xa, Xb, std=TRUE, bootstrap=FALSE)
  # apply(res2,2, quantile, c(0.025, 0.5, 0.975))
  
  round(data.frame(sept_oct=res1[1,], oct_nov=res2[1,]), 2)
  #            sept_oct oct_nov
  # adj.plates     0.29   -0.21
  # eda_dom       -0.26    0.04
  
  # Rennison et al results
  # .34     -.21
  # -.42     .26

  mv_selection_differential <- function(Xa, Xb, std=TRUE, bootstrap=FALSE, niter=4999){
    # Xa contains the data before selection
    # Xb contains the data after selection
    inc_A <- 1:nrow(Xa)
    inc_B <- 1:nrow(Xb)
    if(bootstrap==FALSE){niter <- 0}
    B <- matrix(NA, nrow=(niter+1), ncol=ncol(Xa)) # B is the matrix of b vectors for each iteration of the bootstrap
    colnames(B) <- colnames(Xa)
    for(iter in 1:(niter+1)){
      S <- cov(Xa[inc_A,])
      xbar <- apply(Xb[inc_B], 2, mean) - apply(Xa[inc_A,], 2, mean)
      B[iter,] <- solve(S)%*%xbar
      sd <- sqrt(diag(S))
      if(std==TRUE){B[iter,] <- B[iter,]*sd}
      inc_A <- sample(1:nrow(Xa), replace=TRUE)
      inc_B <- sample(1:nrow(Xb), replace=TRUE)
    }
    return(B)	
  }
  