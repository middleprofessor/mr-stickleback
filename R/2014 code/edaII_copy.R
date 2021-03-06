# reboot Aug 9 2013
library(corpcor) # cov.shrink
library(data.table)

repeatrow <- function(r,n){
	#r is length p
	#creates a n x p matrix with all rows = r
	p <- length(r)
	m <- t(matrix(rep(r,n),ncol=n,nrow=p))
	return(m)
}

repeatcol <- function(c,p){
	#c is length n
	#creates a n x p matrix with all cols = c
	n <- length(c)
	m <- matrix(rep(c,p),ncol=p,nrow=n)
	return(m)
}

center.by.control.means <- function(y,g,tlevel){
	# the data are comprised of multiple fish from 6 families and the
	# families were split into treatment and control levels.
	# to remove family effects, need to center data by family means but
	# use the family mean within control level as this is the best
	# estimate of the mean without selection
	
	# y is the matrix of variables to center (in columns)
	# g is the family (with rows matching that in y)
	# tlevel is the treatment level (1=treated, 0=control)

	ngroups <- length(levels(g)) #number of groups
	n <- dim(y)[1]
	p <- dim(y)[2]
	inc <- which(tlevel==0)
	x <- data.table(Y=y[inc,],family=g[inc])
	xbar <- data.frame(x[,lapply(.SD,mean,na.rm=TRUE),by=family])[,2:(p+1)] #group means
	if(p==1){xbar <- matrix(xbar,ncol=1)}
	yc <- data.frame(matrix(0,nrow=n,ncol=p))
	for(i in 1:n){
		yc[i,] <- y[i,] - xbar[as.numeric(g[i]),]
	}
	#test
#	yc <- data.table(yc,g=g)
#	test <- yc[,lapply(.SD,mean,na.rm=TRUE),by=g] #group means

	return(yc)
}

centerbygroup <- function(y,g){
	#y is a matrix
	#g is grouping factor
	ngroups <- length(levels(g)) #number of groups
	n <- dim(y)[1]
	p <- dim(y)[2]
	x <- data.table(y,g)
	xbar <- data.frame(x[,lapply(.SD,mean,na.rm=TRUE),by=g])[,2:(p+1)] #group means
	if(p==1){xbar <- matrix(xbar,ncol=1)}
	yc <- data.frame(matrix(0,nrow=n,ncol=p))
	for(i in 1:n){
		yc[i,] <- y[i,] - xbar[as.numeric(g[i]),]
	}
	#test
#	yc <- data.table(yc,g=g)
#	test <- yc[,lapply(.SD,mean,na.rm=TRUE),by=g] #group means

	return(yc)
}


return_residuals <- function(y,x){
	# return residuals of y on x
	res <- resid(lm(y~x,na.action=na.exclude))
	return(res)
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

mendelian_randomization <- function(x,y,z,bootstrap=FALSE,niter=4999,A=1,B=2){
	# computes MR of a selection differential, B - A, Not regression!
	# B should be the selected group, A should be the control
	# x is a continuous morphological trait
	# y is NOT survival but experimental treatment
	# z is genotype score
	# A and B are group identification. The differential is B - A
	# if bootstrap then return bootstrap distribution of regression
	N <- length(x)
	groups <- levels(y)
	A_cases <- which(y==groups[A]) # control group cases
	B_cases <- which(y==groups[B]) # selected group cases
	beta.xy <- mean(x[B_cases],na.rm=TRUE) - mean(x[A_cases],na.rm=TRUE)
	beta.zy <- mean(z[B_cases],na.rm=TRUE) - mean(z[A_cases],na.rm=TRUE)
	beta.zx <- summary(lm(x[A_cases]~z[A_cases]))$coefficients[2]
	mr <- beta.zy/beta.zx
	out <- data.frame(beta.xy=beta.xy,beta.zy=beta.zy,beta.zx=beta.zx,mr=mr)
	
	if(bootstrap==TRUE){
		for(i in 1:niter){
			inc <- sample(1:N,N,replace=TRUE)
			x.r <- scale(x[inc])
			y.r <- y[inc]
			z.r <- scale(z[inc])
			A_cases.r <- which(y.r==groups[A]) # control group cases
			B_cases.r <- which(y.r==groups[B]) # selected group cases
			beta.xy <- mean(x.r[B_cases.r],na.rm=TRUE) - mean(x.r[A_cases.r],na.rm=TRUE)
			beta.zy <- mean(z.r[B_cases.r],na.rm=TRUE) - mean(z.r[A_cases.r],na.rm=TRUE)
			beta.zx <- summary(lm(x.r[A_cases.r]~z.r[A_cases.r]))$coefficients[2]
			mr <- beta.zy/beta.zx
			out <- rbind(out,data.frame(beta.xy=beta.xy,beta.zy=beta.zy,beta.zx=beta.zx,mr=mr))
		}
	}
	return(out)
}

TSLS <- function(x,y,z,bootstrap=FALSE,niter=4999,A=1,B=2){
	# computes MR of a selection differential, B - A, Not regression!
	# B should be the selected group, A should be the control
	# x is a continuous morphological trait
	# y is NOT survival but experimental treatment
	# z is genotype score
	# A and B are group identification. The differential is B - A
	# if bootstrap then return bootstrap distribution of regression
	N <- length(x)
	groups <- levels(y)
	A_cases <- which(y==groups[A]) # control group cases
	B_cases <- which(y==groups[B]) # selected group cases
	beta.zx <- summary(lm(x[A_cases]~z[A_cases]))$coefficients[1:2,1]
	xhat <- beta.zx[1] + beta.zx[2]*z
	sx <- sd(xhat[A_cases],na.rm=TRUE)
	sz <- sd(z[A_cases],na.rm=TRUE)
	mr <- (mean(xhat[B_cases],na.rm=TRUE) - mean(xhat[A_cases],na.rm=TRUE))
	beta.zy <- (mean(z[B_cases],na.rm=TRUE) - mean(z[A_cases],na.rm=TRUE))/sz
	return(mr)
}

twotail <- function(x){
	n <- length(x)
	o <- length(which(abs(x) >= abs(x[1])))
	P <- o/n
	return(P)
}

myttest <- function(x,g){
	# paired t test of values in x grouped by g
	a <- x[g==1]
	b <- x[g==0]
	P <- t.test(a,b,paired=TRUE,na.rm=TRUE)$p.value
	return(P)
}

table1 <- function(mydata){
	# computes the table of
	# 1) correlations between eda and trait
	# 2) standardized selection differential on the trait
	# 3) standaridized beta on the trait 
	# using among-individual data centered by family-means
	mydata <- data.table(read.table('eda.txt',header=TRUE,sep='\t'))	
	niter <- 10
	t1 <- array(0,dim=c(niter,5,3))
	bzy <- numeric(niter)
	include.all <- FALSE # if true then include both treatments in computation of correlations but not std dev of X or Z alone
	samples <- 1:dim(mydata)[1]
	for(i in 1:niter){
		data <- mydata[samples,list(family,standlength,ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]
		data <- data[,treat_level:=sapply(mydata[samples,treatment], function(x) if(x=='pred'){1}else{0})]
		data <- data[,eda.score:=as.numeric(mydata[samples,eda.genotype])]

		treated_cases <- which(data[,treat_level]==1)
		control_cases <- which(data[,treat_level]==0)
		# center by family using treatment=0 (control) only
		eda.scores.centered <- center.by.control.means(matrix(data[,eda.score],ncol=1),data[,family],data[,treat_level])[,1] #centered
		#eda.scores.centered <- centerbygroup(matrix(data[,eda.score],ncol=1),data[,family])[,1] #centered
		# standardize by control group
		eda.scores.std <- eda.scores.centered/sd(eda.scores.centered[control_cases],na.rm=TRUE)
		
		# complex includes
			# ant.dorspine
			# sec.dorspine
			# pelspine.len
			# pelgirdle.len
			# plate.number
		# center complex by family means
		sl.centered <- center.by.control.means(as.matrix(data[,standlength]),data[,family],data[,treat_level])[,1]
		complex.centered <- center.by.control.means(as.matrix(data[,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]),data[,family],data[,treat_level])
		# remove covariation due to size using residuals from regression on size
		complex.size_free <- sapply(complex.centered,return_residuals,x=sl.centered)
		# standardize by control sd
		complex.size_free.std <- complex.size_free/repeatrow(apply(complex.size_free[control_cases,],2,sd),dim(complex.size_free)[1])
		colnames(complex.size_free.std) <- c('ant.dorspine','sec.dorspine','pelspine.len','pelgirdle.len','plate.number')
		
		# combine eda and traits in complex
		std_data <- cbind(eda=eda.scores.std, complex.size_free.std)
		
		inc <- control_cases
		if(include.all==TRUE){inc <- c(treated_cases,control_cases)}
		# correlation of eda with individual traits in complex
		which_traits <- c(2,3,4,5,6)
		r <- cor(std_data[inc,],use='pairwise.complete.obs')
		cor_eda_traits <- r[which_traits,1]
		
		# compute vector of standardized selection differentials
		which_traits <- c(2,3,4,5,6)
		dxy <- apply(std_data,2,selection_differential,y=data[,treat_level])[which_traits]

#		xtx <- cov(std_data[inc,which_traits])
		xtx <- cov.shrink(std_data[inc,which_traits],verbose=FALSE)
		xtxi <- solve(xtx)
		bxy <- xtxi%*%dxy
		
		bzy[i] <- (mean(std_data[treated_cases,'eda'],na.rm=TRUE) - mean(std_data[control_cases,'eda'],na.rm=TRUE))

		t1[i,,1] <- cor_eda_traits
		t1[i,,2] <- dxy
		t1[i,,3] <- bxy
		sample1 <- sample(which(mydata[,treatment]=='no'),replace=TRUE)
		sample2 <- sample(which(mydata[,treatment]=='pred'),replace=TRUE)
		samples <- c(sample1,sample2) 
	}
	dimnames(t1) <- list(NULL,colnames(std_data)[which_traits],c('cor_eda','dxy','bxy'))
	meds <- apply(t1,c(2,3),median)
	emp <- t1[1,,]
	low95 <- apply(t1,c(2,3),quantile,prob=0.025)
	up95 <- apply(t1,c(2,3),quantile,prob=0.975)
	low90 <- apply(t1,c(2,3),quantile,prob=0.05)
	up90 <- apply(t1,c(2,3),quantile,prob=0.95)
	zest <- bzy[1]
	zdis <- quantile(bzy,prob=c(0.025,0.5,0.975))
}

table2 <- function(mydata){
	# computes the table of
	# 1) correlations between eda and trait
	# 2) standardized selection differential on the trait
	# 3) standaridized beta on the trait 
	# within each family (so table 1 done 6 times)
	bootit <- FALSE # obtain bootstrap intervals
	permuteit <- TRUE # obtain permuted p-values
	mydata <- mydata[,family:=factor(family)]
	niter <- 1000
	t1 <- array(0,dim=c(niter,5,4))
	bzy <- numeric(niter)
	include.all <- FALSE # if true then include both treatments in computation of correlations but not std dev of X or Z alone
	samples <- 1:dim(mydata)[1]
	for(i in 1:niter){
		data <- mydata[samples,list(family,standlength,ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]
		data <- data[,treat_level:=sapply(mydata[samples,treatment], function(x) if(x=='pred'){1}else{0})]
		data <- data[,eda.score:=as.numeric(mydata[samples,eda.genotype])]
		if(permuteit==TRUE & i>1){
			data[,treat_level:=sample(data[,treat_level])]
		}

		treated_cases <- which(data[,treat_level]==1)
		control_cases <- which(data[,treat_level]==0)
		
		# center by family using treatment=0 (control) only
		eda.scores.centered <- center.by.control.means(matrix(data[,eda.score],ncol=1),data[,family],data[,treat_level])[,1] #centered
		#eda.scores.centered <- centerbygroup(matrix(data[,eda.score],ncol=1),data[,family])[,1] #centered
		# standardize by control group
		eda.scores.std <- eda.scores.centered/sd(eda.scores.centered[control_cases],na.rm=TRUE)

		# complex includes
			# ant.dorspine
			# sec.dorspine
			# pelspine.len
			# pelgirdle.len
			# plate.number
		# center complex by family means
		sl.centered <- center.by.control.means(as.matrix(data[,standlength]),data[,family],data[,treat_level])[,1]
		complex.centered <- center.by.control.means(as.matrix(data[,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]),data[,family],data[,treat_level])
		# remove covariation due to size using residuals from regression on size
		complex.size_free <- sapply(complex.centered,return_residuals,x=sl.centered)
		# standardize by control sd
		complex.size_free.std <- complex.size_free/repeatrow(apply(complex.size_free[control_cases,],2,sd),dim(complex.size_free)[1])
		colnames(complex.size_free.std) <- c('ant.dorspine','sec.dorspine','pelspine.len','pelgirdle.len','plate.number')
		
		# combine eda and traits in complex
		std_data <- cbind(eda=eda.scores.std, complex.size_free.std)
		
		inc <- control_cases
		if(include.all==TRUE){inc <- c(treated_cases,control_cases)}
		# correlation of eda with individual traits in complex
		which_traits <- c(2,3,4,5,6)
		r <- cor(std_data[inc,],use='pairwise.complete.obs')
		cor_eda_traits <- r[which_traits,1]
		
		# compute vector of standardized selection differentials
		which_traits <- c(2,3,4,5,6)
		dxy.by.fam <- matrix(0,nrow=length(levels(data[,family])),ncol=length(which_traits))
		colnames(dxy.by.fam) <- colnames(std_data)[which_traits]
		for(which_family in 1:length(levels(data[,family]))){
			rows <- which(data[,family]==levels(data[,family])[which_family])
			dxy.by.fam[which_family,] <- apply(std_data[rows,],2,selection_differential,y=data[rows,treat_level],std=TRUE)[which_traits]
		}
		dxy.med <- apply(dxy.by.fam,2,median) #med among families
		dxy.bar <- apply(dxy.by.fam,2,mean) #mean among families
		dxy <- apply(std_data,2,selection_differential,y=data[,treat_level])[which_traits] #among individuals that are mean centered by families


#		xtx <- cov(std_data[inc,which_traits])
		xtx <- cov.shrink(std_data[inc,which_traits],verbose=FALSE)
		xtxi <- solve(xtx)
		bxy <- xtxi%*%dxy
		
		bzy[i] <- (mean(std_data[treated_cases,'eda'],na.rm=TRUE) - mean(std_data[control_cases,'eda'],na.rm=TRUE))

		t1[i,,1] <- cor_eda_traits
		t1[i,,2] <- dxy
		t1[i,,3] <- dxy.med
		t1[i,,4] <- dxy.bar
		if(bootit==TRUE){			
			sample1 <- sample(which(mydata[,treatment]=='no'),replace=TRUE)
			sample2 <- sample(which(mydata[,treatment]=='pred'),replace=TRUE)
			samples <- c(sample1,sample2) # this will give bootstrap estimates
		}
	}
	dimnames(t1) <- list(NULL,colnames(std_data)[which_traits],c('cor_eda','dxy-ind','dxy-med','dxy-bar'))
	meds <- apply(t1,c(2,3),median)
	emp <- t1[1,,]
	low95 <- apply(t1,c(2,3),quantile,prob=0.025,na.rm=TRUE)
	up95 <- apply(t1,c(2,3),quantile,prob=0.975,na.rm=TRUE)
	low90 <- apply(t1,c(2,3),quantile,prob=0.05)
	up90 <- apply(t1,c(2,3),quantile,prob=0.95)
	zest <- bzy[1]
	zdis <- quantile(bzy,prob=c(0.025,0.5,0.975))
}

table3 <- function(mydata){
	# attempts to recover table 1 from Marchinko
	# 1) correlations between eda and trait
	# 2) standardized selection differential on the trait
	# 3) standaridized beta on the trait 
	# within each family (so table 1 done 6 times)
	bootit <- TRUE # obtain bootstrap intervals
	permuteit <- FALSE # obtain permuted p-values
	scale_array <- c(TRUE,FALSE) #scale within family or ignoring family
	zeroes_array <- c(TRUE,FALSE) # analyze zeros or convert to NA
	ntests <- length(scale_array)*length(zeroes_array)
	test <- c(0,0)
	testarray <- t(matrix(c(c(1,1),c(1,2),c(2,1),c(2,2)),ncol=4))
	table3 <- NULL
	table3.p <- NULL
	table3.lo95 <- NULL
	table3.hi95 <- NULL
	for(which_test in 1:ntests){
		mydata <- data.table(read.table('eda.txt',header=TRUE,sep='\t'))
		mydata <- mydata[,family:=factor(family)]
		test <- testarray[which_test,]
		scale_by_family <- scale_array[test[1]]
		include_zeroes <- zeroes_array[test[2]]
		#replace zeds
		if(include_zeroes==FALSE){
			tdata <- mydata[,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len)]
			tdata[tdata==0] <- NA
			mydata <- mydata[,ant.dorspine:=tdata[,ant.dorspine]]
			mydata <- mydata[, sec.dorspine:=tdata[, sec.dorspine]]
			mydata <- mydata[, pelspine.len:=tdata[, pelspine.len]]
			mydata <- mydata[, pelgirdle.len:=tdata[, pelgirdle.len]]
		}
		
		niter <- 1000
		t1 <- array(0,dim=c(niter,5,4))
		bzy <- numeric(niter)
		include.all <- FALSE # if true then include both treatments in computation of correlations but not std dev of X or Z alone
		samples <- 1:dim(mydata)[1]
	
		#declare matrices
		armor_traits <- c('ant.dorspine','sec.dorspine','pelspine.len','pelgirdle.len','plate.number')
		p <- length(armor_traits)
		dxy.by.fam <- matrix(0,nrow=length(levels(data[,family])),ncol=p)
		treated.bar <- matrix(0,nrow=length(levels(data[,family])),ncol=p)
		treated.med <- matrix(0,nrow=length(levels(data[,family])),ncol=p)
		untreated.bar <- matrix(0,nrow=length(levels(data[,family])),ncol=p)
		untreated.med <- matrix(0,nrow=length(levels(data[,family])),ncol=p)
		colnames(dxy.by.fam) <- armor_traits
		colnames(treated.bar) <- armor_traits
		colnames(treated.med) <- armor_traits
		colnames(untreated.bar) <- armor_traits
		colnames(untreated.med) <- armor_traits
	
	
		for(i in 1:niter){
			data <- mydata[samples,list(family,standlength,ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]
			data <- data[,treat_level:=sapply(mydata[samples,treatment], function(x) if(x=='pred'){1}else{0})]
			data <- data[,eda.score:=as.numeric(mydata[samples,eda.genotype])]
			
			# permute if i > 1
			if(permuteit==TRUE & i>1){
				data[,treat_level:=sample(data[,treat_level])]
			}
	
			treated_cases <- which(data[,treat_level]==1)
			control_cases <- which(data[,treat_level]==0)
			
			#scale by all individuals not within families
			armor.all <-sapply(data[,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)],return_residuals,x=data[,standlength])
			
			#center by family (after scaling) and get selection differential
			if(scale_by_family==FALSE){
				armor.all.c <- center.by.control.means(as.matrix(armor.all),data[,family],data[,treat_level])		
				dxy <- apply(armor.all.c,2,selection_differential,y=data[,treat_level],std=TRUE)
				#for bxy need to de-standardize dxy using dxy*sdx and then re-standardize
				sdx <- apply(armor.all.c[data[,treat_level]==0,],2,sd,na.rm=TRUE)
				bxy <- solve(cov(armor.all.c,use='pairwise.complete.obs'))%*%(dxy*sdx)*sdx

			}
			
			#center by family (including SL) and get selection differential
			if(scale_by_family==TRUE){
				sl.c <- center.by.control.means(as.matrix(data[,standlength]),data[,family],data[,treat_level])[,1]
				armor.all.c <- center.by.control.means(as.matrix(data[,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]),data[,family],data[,treat_level])
				# remove covariation due to size using residuals from regression on size
				armor.all.c <- sapply(armor.all.c,return_residuals,x=sl.c)
				# get selection differential standardized by control sd
				dxy <- apply(armor.all.c,2,selection_differential,y=data[,treat_level],std=TRUE)
				sdx <- apply(armor.all.c[data[,treat_level]==0,],2,sd,na.rm=TRUE)
				bxy <- solve(cov(armor.all.c,use='pairwise.complete.obs'))%*%(dxy*sdx)*sdx

			}
			
			for(which_family in 1:length(levels(data[,family]))){
				rows <- which(data[,family]==levels(data[,family])[which_family])
				
				if(scale_by_family==TRUE){
					armor <-sapply(data[rows,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)],return_residuals,x=data[rows,standlength])
				}else{armor <- armor.all[rows,]}
	
				dxy.by.fam[which_family,] <- apply(armor,2,selection_differential,y=data[rows,treat_level],std=TRUE)
				
				if(i==1){
					treated.bar[which_family,] <- apply(armor[data[rows,treat_level]==1,],2,mean,na.rm=TRUE)
					treated.med[which_family,] <- apply(armor[data[rows,treat_level]==1,],2,median,na.rm=TRUE)
					untreated.bar[which_family,] <- apply(armor[data[rows,treat_level]==0,],2,mean,na.rm=TRUE)
					untreated.med[which_family,] <- apply(armor[data[rows,treat_level]==0,],2,median,na.rm=TRUE)			
				}
			}
	
			dxy.med <- apply(dxy.by.fam,2,median) #med among families
			dxy.bar <- apply(dxy.by.fam,2,mean) #mean among families
			
			t1[i,,1] <- bxy
			t1[i,,2] <- dxy
			t1[i,,3] <- dxy.med
			t1[i,,4] <- dxy.bar
			if(bootit==TRUE){			
				sample1 <- sample(which(mydata[,treatment]=='no'),replace=TRUE)
				sample2 <- sample(which(mydata[,treatment]=='pred'),replace=TRUE)
				samples <- c(sample1,sample2) # this will give bootstrap estimates
			}
		}
		
		dimnames(t1) <- list(NULL,colnames(std_data)[which_traits],c('bxy','dxy-ind','dxy-med','dxy-bar'))
		emp <- t1[1,,]

		p.ttest.bar <- apply(rbind(treated.bar,untreated.bar),2,myttest,g=c(rep(1,6),rep(0,6)))
		p.ttest.med <- apply(rbind(treated.med,untreated.med),2,myttest,g=c(rep(1,6),rep(0,6)))
	
		if(permuteit==TRUE){
			p.permute <- apply(t1,c(2,3),twotail)
			rownames(p.permute) <- rownames(emp)
			colnames(p.permute) <- colnames(emp)
		}
		
		if(bootit==TRUE){			
			lo95 <- apply(t1,c(2,3),quantile,prob=0.025,na.rm=TRUE)
			hi95 <- apply(t1,c(2,3),quantile,prob=0.975,na.rm=TRUE)
#			low90 <- apply(t1,c(2,3),quantile,prob=0.05,na.rm=TRUE)
#			up90 <- apply(t1,c(2,3),quantile,prob=0.95,na.rm=TRUE)
		}

		table3 <- rbind(table3,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='bxy',t(t1[1,,1])))
		table3 <- rbind(table3,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy',t(t1[1,,2])))
		table3 <- rbind(table3,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy-bar',t(emp[,4])))
		table3 <- rbind(table3,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy-med',t(emp[,3])))

		if(permuteit==TRUE){
			table3.p <- rbind(table3.p,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='p.permute.bxy',t(p.permute[,1])))
			table3.p <- rbind(table3.p,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='p.permute.dxy',t(p.permute[,2])))
			table3.p <- rbind(table3.p,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='p.ttest.bar',t(p.ttest.bar)))
			table3.p <- rbind(table3.p,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='p.ttest.med',t(p.ttest.med)))
			table3.p <- rbind(table3.p,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='p.permute.bar',t(p.permute[,4])))
			table3.p <- rbind(table3.p,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='p.permute.med',t(p.permute[,3])))
		}
		if(bootit==TRUE){
			table3.hi95 <- rbind(table3.hi95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='bxy', t(hi95[,1])))
			table3.hi95 <- rbind(table3.hi95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy', t(hi95[,2])))
			table3.hi95 <- rbind(table3.hi95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy-bar', t(hi95[,3])))
			table3.hi95 <- rbind(table3.hi95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy-med', t(hi95[,4])))

			table3.lo95 <- rbind(table3.lo95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='bxy',t(lo95[,1])))
			table3.lo95 <- rbind(table3.lo95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy',t(lo95[,2])))
			table3.lo95 <- rbind(table3.lo95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy-bar', t(lo95[,3])))
			table3.lo95 <- rbind(table3.lo95,data.frame(scale_by_family=scale_by_family,zeroes=include_zeroes,test='dxy-med', t(lo95[,4])))		
		}

	}
	
	if(permuteit==TRUE){	
		write.table(table3,'table3.permute.txt',quote=FALSE,sep='\t',row.names=FALSE)
		write.table(table3.p,'table3.p-values.txt',quote=FALSE,sep='\t',row.names=FALSE)
	}
	if(bootit==TRUE){	
		write.table(table3,'table3.boot.txt',quote=FALSE,sep='\t',row.names=FALSE)
		write.table(table3.hi95,'table3.hi95.txt',quote=FALSE,sep='\t',row.names=FALSE)
		write.table(table3.lo95,'table3.lo95.txt',quote=FALSE,sep='\t',row.names=FALSE)
	}

}

doit <- function(){
	
	#Start Here!
	mydata <- data.table(read.table('eda.txt',header=TRUE,sep='\t'))
	# columns are
		# population
		# family
		# treatment {no = -insect, pred = +insect}
		# individual
		# standlength
		# ant.dorspine
		# sec.dorspine
		# pelspine.len
		# pelgirdle.len
		# plate.number
		# eda.genotype {AA, Aa, aa}
		# pelgirdle.presence {1=yes, 0 = no}
		# pelspine.presence {1 = yes, 0 = no}
	
	# get control and treated cases
	treated_cases <- which(mydata[,treatment]=='pred')
	control_cases <- which(mydata[,treatment]=='no')

	
	# transform genotype into score
	mydata <- mydata[,eda.score:=as.numeric(eda.genotype)]
	# center by family
	eda.scores.centered <- centerbygroup(matrix(mydata[,eda.score],ncol=1),mydata[,family])[,1] #centered
	# standardize by control group
	eda.scores.std <- eda.scores.centered/sd(eda.scores.centered[control_cases],na.rm=TRUE)
	
	# complex includes
		# ant.dorspine
		# sec.dorspine
		# pelspine.len
		# pelgirdle.len
		# plate.number
	# center complex by family means
	sl.centered <- centerbygroup(as.matrix(mydata[,standlength]),mydata[,family])[,1]
	complex.centered <- centerbygroup(as.matrix(mydata[,list(ant.dorspine,sec.dorspine,pelspine.len,pelgirdle.len,plate.number)]),mydata[,family])
	# remove covariation due to size using residuals from regression on size
	complex.size_free <- sapply(complex.centered,return_residuals,x=sl.centered)
	# standardize by control sd
	complex.size_free.std <- complex.size_free/repeatrow(apply(complex.size_free[control_cases,],2,sd),dim(complex.size_free)[1])
	colnames(complex.size_free.std) <- c('ant.dorspine','sec.dorspine','pelspine.len','pelgirdle.len','plate.number')
	
	# combine eda and traits in complex
	std_data <- cbind(eda=eda.scores.std, complex.size_free.std)
	
	# correlation of eda with individual traits in complex
	which_rows <- c(control_cases,treated_cases)
	which_traits <- c(2,3,4,5,6)
	r <- cor(std_data[which_rows,],use='pairwise.complete.obs')
	cor_eda_traits <- r[which_traits,1]
	
	# compute vector of standardized selection differentials
	which_traits <- c(2,3,4,5,6)
	dxy <- apply(std_data,2,selection_differential,y=mydata[,treatment])[which_traits]
	n <- length(control_cases)
	xtx <- t(std_data[control_cases,which_traits])%*%std_data[control_cases,which_traits]/(n-1) #not = to cor_eda_traits because that threw out incomplete rows
	xtxi <- solve(xtx)
	bxy <- xtxi%*%dxy
	bxy.mr <- do.call(rbind,apply(complex.size_free.std,2,mendelian_randomization,y=mydata[,treatment], z=eda.scores.std))
	bxy.tsls <- apply(complex.size_free.std,2,TSLS,y=mydata[,treatment],z=eda.scores.std)
	table1 <- data.frame(cor_eda=cor_eda_traits,dxy=dxy,bxy=bxy,bxy.mr=bxy.mr[,'mr'])
	
	#bootstrap table 1

	# compute composite armor score as the weighted sum of the complex
	weights.1 <- c(1,1,1,1,1) # all traits 
	# not that plates and dorsal spines are correlated with each other and with eda but pelvic
	# complex much less so, this suggests a weighting is...
	weights.2 <- c(1,1,0,0,1) # dorsal spines & lateral plates but not pelvic complex
	# weight as PC1
	inc <- c(1,2,5)
	pc <- eigen(cor(complex.size_free.std[,inc]))$vectors[1:3,1]
	weights.3 <- weights.2
	weights.3[1:2] <- pc[1:2]
	weights.3[5] <- pc[3]
	
	armor_score.2 <- complex.size_free.std%*%weights.2
	armor_score.2 <- armor_score.2/sd(armor_score.2[control_cases]) 
	cor_eda_composite.2 <- cor(eda.scores.std[control_cases],armor_score.2[control_cases],use='complete.obs')
	Fstat <- summary(lm(armor_score.2[control_cases]~eda.scores.std[control_cases]))$fstatistic
	beta.xy.mr.2 <- mendelian_randomization(armor_score.2,mydata[,treatment],eda.scores.std,TRUE,niter=4999)
	beta.xy.mr.2.quants <- sapply(beta.xy.mr.2,quantile,c(0.025,0.975))
	cor_composite_traits <- cor(data.frame(armor=armor_score.2, complex.size_free.std)[control_cases,],use='complete.obs')

	armor_score.3 <- complex.size_free.std%*%weights.3
	armor_score.3 <- armor_score.3/sd(armor_score.3[control_cases]) 
	cor_eda_composite.3 <- cor(eda.scores.std[control_cases],armor_score.3[control_cases],use='complete.obs')
	Fstat <- summary(lm(armor_score.3[control_cases]~eda.scores.std[control_cases]))$fstatistic
	beta.xy.mr.3 <- mendelian_randomization(armor_score.3,mydata[,treatment],eda.scores.std,TRUE,niter=4999)
	beta.xy.mr.3.quants <- sapply(beta.xy.mr.3,quantile,c(0.025,0.975))
	cor_composite_traits <- cor(data.frame(armor=armor_score.3, complex.size_free.std)[control_cases,],use='complete.obs')

	#individual parts of MR
	beta.xy <- mean(armor_score.2[treated_cases],na.rm=TRUE) - mean(armor_score.2[control_cases],na.rm=TRUE)
	beta.zy <- mean(eda.scores.std[treated_cases],na.rm=TRUE) - mean(eda.scores.std[control_cases],na.rm=TRUE)
	beta.zx <- summary(lm(armor_score.2[control_cases]~eda.scores.std[control_cases]))$coefficients[2]
	cor.zy <- cor(armor_score.2[control_cases],eda.scores.std[control_cases],use='complete.obs')
	mr <- beta.zy/beta.zx
	
	
	armor_score.2 <- complex.size_free.std%*%weights.2
	armor_score.2 <- armor_score.2/sd(armor_score.2[control_cases]) 
	cor_eda_composite.2 <- cor(eda.scores.std[control_cases], armor_score.2[control_cases],use='complete.obs')
	beta.xy.mr.2 <- mendelian_randomization(armor_score.2,mydata[,treatment],eda.scores.std,FALSE)

	armor_score.3 <- complex.size_free.std%*%weights.3
	armor_score.3 <- armor_score.3/sd(armor_score.3[control_cases]) 
	cor_eda_composite.3 <- cor(eda.scores.std[control_cases],armor_score.3[control_cases],use='complete.obs')
	beta.xy.mr.3 <- mendelian_randomization(armor_score.3,mydata[,treatment],eda.scores.std,FALSE)

}

testpcalg <- function(){
	data <- cbind(std_data,treat=as.numeric(mydata[,treatment]))
	data <- na.omit(data)
	suffStat <- list(C=cor(data),n=nrow(data))
	pc.fit <- pc(suffStat,indepTest=gaussCItest,p=ncol(data),alpha=0.05)
	plot(pc.fit,main="")
}