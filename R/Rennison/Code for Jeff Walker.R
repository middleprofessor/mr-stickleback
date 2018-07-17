#Size adjustment of plates by breakpoint regression. 

# Preparation steps:
# 1. breakpoint regression of log plates using fish from Sept, Oct, Nov, Jan (breakpoint is 34.0)
# 2. Calculation of residuals
# 3. Give everyone a length of 34mm and predict their logplates by their genotype
# 4. Add the residuals back to the prediction to yield adusted logplates. 
# 5. Take exp(adj log plates) to yield adjusted plates
# 6. Set plates to 32 if plates > 32. 
#
# as follows:
# ----------
x <- read.csv("data.csv", stringsAsFactors=FALSE)


# Segmented regression logplates without individual ponds
# log plates, full model: 3 genotypes, separate slopes, separate intercepts
z <- nls(logplates ~  (a + b * length) * as.numeric(genotype == "CC" & length<=k) + 
					  (a + b * k) * as.numeric(genotype == "CC" & length>k) +
					  (c + d * length) * as.numeric(genotype == "CL" & length<=k) +
					  (c + d * k) * as.numeric(genotype == "CL" & length>k) +
					  (e + f * length) * as.numeric(genotype == "LL" & length<=k) +
					  (e + f * k) * as.numeric(genotype == "LL" & length>k),
					  data=x, algorithm = "default",
					  start=list(a=1, c=1, e=1, b=0.1, d=0.1, f=0.1, k=34), control=list(warnOnly=TRUE))

# Calculate and save residuals
res <- residuals(z)
# check that length of res equals the number of rows of x
length(res)
nrow(x)

# Calculate an adjusted plate number for 34mm, constrained to a maximum of 32
x1 <- cbind.data.frame(genotype = c("CC","CL","LL"), length=c(34,34,34), stringsAsFactors=FALSE)
z1 <- predict(z, newdata=x1)
adjlogplates <- z1[match(x$genotype, c("CC","CL","LL"))]
adjlogplates <- adjlogplates + res
adjplates <- exp(adjlogplates)
# cap number of plates at 32
adjplates[adjplates > 32] <- 32
x$adjplates <- adjplates


### Lande-Arnold analysis for estimating partial selection coefficients
x <- read.csv("data.csv", stringsAsFactors=FALSE)

#Function for generating selection coefficients in the lande arnold style 
beta.before.after <- function(means.before, means.after, var.before){
	# means.before, means.after are vectors of the same length
	# var.before is a symmetric covariance matrix with dimensions equal to length of mean vectors
	#
	if(!all(c(dim(var.before), length(means.after)) == length(means.before)))
		stop("dimensions of vectors and matrix do not agree")
	means.after <- unlist(as.vector(means.after))
	means.before <- unlist(as.vector(means.before))
	results <- list()
	results$s <- means.after - means.before
	results$beta <- solve(var.before) %*% results$s
	colnames(results$beta) <- "beta"
	results$beta.std <- results$beta*sqrt(diag(var.before))
	return(results)
	}

#Two genotype analysis of Lateral plate phenotype and Eda genotype

#Subsetting genotypes
CC <- subset(x,x$genotype=="CC")
CL <- subset(x,x$genotype=="CL")

#Subsetting Genotype by month and inspecting the means
CCs <- subset(CC,CC$month=="Sep-06")
#mean(CCs$adj.plates)
#22.1

CLs <- subset(CL,CL$month=="Sep-06")
#mean(CLs$adj.plates)
#19.74286

CCo <- subset(CC,CC$month=="Oct-06")
#mean(CCo$adj.plates) 
#25.82353

CLo <- subset(CL,CL$month=="Oct-06")
#mean(CLo$adj.plates)
#19.14286

CCn <- subset(CC,CC$month=="Nov-06")
#mean(CCn$adj.plates)
#24.31429 

CLn <- subset(CL,CL$month=="Nov-06")
#mean(CLn$adj.plates)
#18.72414


#Defining Genotype Frequncies 

#genotype frequencies from Barrett et al 2008 data
septfreq <- c(0.47, 0.53)
octfreq <- c(0.7, 0.3)
novfreq <- c(0.54, 0.46)


#Define number of individuals to sample 
number <- c(65,55,64)  

#number of repetitions of the loop
n_reps <- 20000


# make 2 matrices to store the selection coeff for sept_oct and oct_nov
oct_nov <- sept_oct <- matrix(NA, nrow=n_reps,ncol=6)
colnames(oct_nov) <- colnames(sept_oct) <-c("s1","s2","beta1","beta2","beta.std1","beta.std2")

#-----


for(i in 1:n_reps) {
	#Sept
	septgeno <- sample (c(0,1), size=200, prob=rep(septfreq), replace=TRUE)
		septgeno2 <- sample(septgeno, size=65, prob=NULL, replace=FALSE )
	CCgeno1 <- subset(septgeno2, septgeno2 ==0)
		sum_CCgeno1 <- length(CCgeno1)
		sadjplates<- vector(length=number[1]) 
		sadjplates[septgeno2==0] <- sample(x = CCs$adj.plates, size = sum_CCgeno1, replace=TRUE, prob = NULL)
		
		
		CLgeno1 <- subset(septgeno2, septgeno2 ==1)
		sum_CLgeno1 <- length(CLgeno1)
		sadjplates[septgeno2==1] <- sample(x = CLs$adj.plates, size = sum_CLgeno1, replace=TRUE, prob = NULL)
	
	
	x1.sept<-cbind(sadjplates, septgeno2) 
	means.sept <- apply(x1.sept, 2, FUN=mean, na.rm=TRUE)
	P.sept <- var(x1.sept)
	
#Oct
	octgeno <- sample (c(0,1), size=200, prob=rep(octfreq), replace=TRUE)
	octgeno2 <- sample(octgeno, size=55, prob=NULL, replace=FALSE )
	CCgeno2 <- subset(octgeno2, octgeno2 ==0)
	sum_CCgeno2 <- length(CCgeno2)
		oadjplates<- vector(length=number[2]) 
		oadjplates[octgeno2==0] <- sample(x = CCo$adj.plates, size = sum_CCgeno2, replace=TRUE, prob = NULL)
		
		
		CLgeno2 <- subset(octgeno2, octgeno2 ==1)
		sum_CLgeno2 <- length(CLgeno2)
		oadjplates[octgeno2==1] <- sample(x = CLo$adj.plates, size = sum_CLgeno2, replace=TRUE, prob = NULL)
		

	x1.oct<-cbind(oadjplates, octgeno2) 
	means.oct <- apply(x1.oct, 2, FUN=mean, na.rm=TRUE)
	P.oct <- var(x1.oct)
#	P.oct[2,2]<- var(octgeno)

sept_oct[i,] <- sept.oct <- unlist(beta.before.after(means.sept, means.oct, P.sept))


	
	#Nov	
		novgeno <- sample (c(0,1), size=200, prob=rep(novfreq), replace=TRUE)
		novgeno2 <- sample(novgeno, size=74, prob=NULL, replace=FALSE )
	CCgeno3 <- subset(novgeno2, novgeno2 ==0)
	sum_CCgeno3 <- length(CCgeno3)
		nadjplates<- vector(length=number[3]) 
		nadjplates[novgeno2==0] <- sample(x = CCn$adj.plates, size = sum_CCgeno3, replace=TRUE, prob = NULL)

		
		
		CLgeno3 <- subset(novgeno2, novgeno2 ==1)
		sum_CLgeno3 <- length(CLgeno3)
		nadjplates[novgeno2==1] <- sample(x = CLn$adj.plates, size = sum_CLgeno3, replace=TRUE, prob = NULL)
		
		
	x1.nov<-cbind(nadjplates, novgeno2) 
	means.nov <- apply(x1.nov, 2, FUN=mean, na.rm=TRUE)
	P.nov <- var(x1.nov)
	#P.nov[2,2]<- var(novgeno)
	oct_nov[i,] <- oct.nov <- unlist(beta.before.after(means.oct, means.nov, P.oct))
	

}

#SeptOct
#Phenotype
mean(sept_oct[,5])
#0.3435019
#Command to retrive the 95% Confidence Interval 
quantile(sept_oct[,5], probs=seq(0,1,0.025))

#Genotype
mean(sept_oct[,6])
#-0.4164707
quantile(sept_oct[,6], probs=seq(0,1,0.025))

#OctNov
#Phenotype
x1[1] <- mean(oct_nov[,5])
#-0.2104259
quantile(oct_nov[,5], probs=seq(0,1,0.025))

#Genotype
mean(oct_nov[,6])
#0.2565881
quantile(oct_nov[,6], probs=seq(0,1,0.025))

#Visualizing the coefficients 
plot(oct_nov[,5] ~ oct_nov[,6])
plot(sept_oct[,5] ~ sept_oct[,6])


#___________________

#Three Genotypes Case

#Subsetting the genotypes 

CC <- subset(x,x$genotype=="CC")
CL <- subset(x,x$genotype=="CL")
LL <- subset(x,x$genotype=="LL")

#Quality control, drop LL fish with plates above 15 
LL <- LL[LL$adj.plates < 15,]

#Subsetting genotypes by month to generate phenotype distributions
CCs <- subset(CC,CC$month=="2006-09")
CCo <- subset(CC,CC$month=="2006-10")
CCn <- subset(CC,CC$month=="2006-11")
CLs <- subset(CL,CL$month=="2006-09")
CLo <- subset(CL,CL$month=="2006-10")
CLn <- subset(CL,CL$month=="2006-11")
LLs <- subset(LL,LL$month=="2006-09")
LLo <- subset(LL,LL$month=="2006-10")
LLn <- subset(LL,LL$month=="2006-11")

#Setting the genotype frequencies 
#revised genotype frequencies from Barrett et al 2008 estimates

septfreq <- c(0.37,0.425,0.205)
octfreq <- c(0.557, 0.233, 0.21)
novfreq <- c(0.41, 0.35, 0.24)


#define number of individuals to sample for each month
number <- c(80,77,84)  

#Setting the number of repetitions of the loop
n_reps <- 10000 

# make 2 matrices to store the selection coeffs for sept_oct and oct_nov
oct_nov <- sept_oct <- matrix(NA, nrow=n_reps,ncol=9)
colnames(oct_nov) <- colnames(sept_oct) <-c("s1","s2","s3","beta1","beta2","beta3","beta.std1","beta.std2","beta.std3")

#The loop for generating the selection coefficients 

for(i in 1:n_reps) {
	#Sept
	septgeno1 <- sample (c(1,0,-1), size=200, prob=rep(septfreq), replace=TRUE)
	septgeno <- sample(septgeno1, size=80, prob=NULL, replace=FALSE )
	CCgeno1 <- subset(septgeno, septgeno ==1)
		sum_CCgeno1 <- length(CCgeno1)
		sadjplates<- vector(length=number[1]) 
		sadjplates[septgeno==1] <- sample(x = CCs$adj.plates, size = sum_CCgeno1, replace=TRUE, prob = NULL)
		
		CLgeno1 <- subset(septgeno, septgeno ==0)
		sum_CLgeno1 <- length(CLgeno1)
		sadjplates[septgeno==0] <- sample(x = CLs$adj.plates, size = sum_CLgeno1, replace=TRUE, prob = NULL)
		
		LLgeno1 <- subset(septgeno, septgeno ==-1)
		sum_LLgeno1 <- length(LLgeno1)
		sadjplates[septgeno==-1] <- sample(x = LLs$adj.plates, size = sum_LLgeno1, replace=TRUE, prob = NULL)
		
sdom <- vector(length=number[1]) 
sdom[septgeno==1] <- 0
sdom[septgeno==0] <- 1
sdom[septgeno==-1] <- 0


	x1.sept<-cbind(sadjplates, septgeno, sdom) 
	means.sept <- apply(x1.sept, 2, FUN=mean, na.rm=TRUE)
	P.sept <- var(x1.sept)

	
#Oct
	octgeno1 <- sample (c(1,0,-1), size=200, prob=rep(octfreq), replace=TRUE)
	octgeno <- sample(octgeno1, size=77, prob=NULL, replace=FALSE )
	CCgeno2 <- subset(octgeno, octgeno ==1)
	sum_CCgeno2 <- length(CCgeno2)
		oadjplates<- vector(length=number[2]) 
		oadjplates[octgeno==1] <- sample(x = CCo$adj.plates, size = sum_CCgeno2, replace=TRUE, prob = NULL)
		
		CLgeno2 <- subset(octgeno, octgeno ==0)
		sum_CLgeno2 <- length(CLgeno2)
		oadjplates[octgeno==0] <- sample(x = CLo$adj.plates, size = sum_CLgeno2, replace=TRUE, prob = NULL)
		
		LLgeno2 <- subset(octgeno, octgeno ==-1)
		sum_LLgeno2 <- length(LLgeno2)
		oadjplates[octgeno==-1] <- sample(x = LLo$adj.plates, size = sum_LLgeno2, replace=TRUE, prob = NULL)
		
odom <- vector(length=number[2]) 
odom[octgeno==1] <- 0
odom[octgeno==0] <- 1
odom[octgeno==-1] <- 0



	x1.oct<-cbind(oadjplates, octgeno, odom) 
	means.oct <- apply(x1.oct, 2, FUN=mean, na.rm=TRUE)
	P.oct <- var(x1.oct)


sept_oct[i,] <- sept.oct <- unlist(beta.before.after(means.sept, means.oct, P.sept))

	
	#Nov	
novgeno1 <- sample (c(1,0,-1), size=200, prob=rep(novfreq), replace=TRUE)
		novgeno <- sample(novgeno1, size=84, prob=NULL, replace=FALSE )
	CCgeno3 <- subset(novgeno, novgeno ==1)
	sum_CCgeno3 <- length(CCgeno3)
		nadjplates<- vector(length=number[3]) 
		nadjplates[novgeno==1] <- sample(x = CCn$adj.plates, size = sum_CCgeno3, replace=TRUE, prob = NULL)
		
		CLgeno3 <- subset(novgeno, novgeno ==0)
		sum_CLgeno3 <- length(CLgeno3)
		nadjplates[novgeno==0] <- sample(x = CLn$adj.plates, size = sum_CLgeno3, replace=TRUE, prob = NULL)
		
		LLgeno3 <- subset(novgeno, novgeno ==-1)
		sum_LLgeno3 <- length(LLgeno3)
		nadjplates[novgeno==-1] <- sample(x = LLn$adj.plates, size = sum_LLgeno3, replace=TRUE, prob = NULL)
		
ndom <- vector(length=number[3]) 
ndom[novgeno==1] <- 0
ndom[novgeno==0] <- 1
ndom[novgeno==-1] <- 0



	x1.nov<-cbind(nadjplates, novgeno, ndom) 
	means.nov <- apply(x1.nov, 2, FUN=mean, na.rm=TRUE)
	P.nov <- var(x1.nov)
#	P.nov[2,2] <- var(novgeno1)
#	P.nov[3,3] <- var(ndom2)
	oct_nov[i,] <- oct.nov <- unlist(beta.before.after(means.oct, means.nov, P.oct))
	

}

#Outputting the results

#SeptOct
#Phenotype
mean(sept_oct[,7])
#0.3478917
quantile(sept_oct[,7], probs=seq(0,1,0.025))
-0.079984794 - 0.784467144

#Genotype
mean(sept_oct[,8])
# -0.03619753
quantile(sept_oct[,8], probs=seq(0,1,0.025))

#Dominance
mean(sept_oct[,9])
# -0.4706104
quantile(sept_oct[,9], probs=seq(0,1,0.025))

#visualizing the additive genotype and phenotype selection estimates
plot(sept_oct[,8] ~ sept_oct[,7])
#visualizing the dominance and phenotype selection estimates
plot(sept_oct[,9] ~ sept_oct[,7])

#OctNov
#Phenotype
mean(oct_nov[,7])
#-0.162349
quantile(oct_nov[,7], probs=seq(0,1,0.025))
-0.830665197 - 0.407847299 

#Genotype
mean(oct_nov[,8])
#-0.03743741
quantile(oct_nov[,8], probs=seq(0,1,0.025))
-0.64617372 - 0.68610132


#Dominance
 mean(oct_nov[,9])
#0.2612135
quantile(oct_nov[,9], probs=seq(0,1,0.025))
-0.136137216 - 0.705515806




#Method for estimating Standardized univarite selection intenesities 

#Univariate estimates of selection on genotype

#Setting the genotype frequncies 

septfreq <- c(0.37,0.425,0.205)
octfreq <- c(0.557, 0.233, 0.21)
novfreq <- c(0.41, 0.35, 0.24)


#Setting the number of repetitions of the loop
n_reps <- 10000

#make 2 matrices to store the selection intensities for sept_oct and oct_nov
oct_nov <- sept_oct <- matrix(NA,nrow=n_reps,ncol=1)
colnames(oct_nov) <- colnames(sept_oct) <-c("intensity")

for(i in 1:n_reps){


	septgeno1 <- sample (c(0,1,0), size=200, prob=rep(septfreq), replace=TRUE)
	octgeno1 <- sample (c(0,1,0), size=200, prob=rep(octfreq), replace=TRUE)
	novgeno1 <- sample (c(0,1,0), size=200, prob=rep(novfreq), replace=TRUE)


	means.sept <- mean(septgeno1)
	sdS <- sd(septgeno1)
	means.oct <- mean(octgeno1)
	sdO <- sd(octgeno1)
	means.nov <- mean(novgeno1)
	sdN <- sd(novgeno1)

		i_SO <- ((means.oct - means.sept)/sdS)
		sept_oct[i,] <- sept.oct <- unlist(i_SO)
		i_ON <- ((means.nov - means.oct)/sdO)
		oct_nov[i,] <- oct.nov <- unlist(i_ON)
}
mean(sept_oct[,1])
#-0.3881128
quantile(sept_oct[,1], probs=seq(0,1,0.025))
#range 2.5% to 97.5%
#-0.55885468 - -0.20844135

mean(oct_nov[,1])
#0.2795994
quantile(oct_nov[,1], probs=seq(0,1,0.025))
0.06629142 - 0.52125003

---------------------

#Univariate estimates of selection on phenotype

#Subsetting of the data 
Sept <- subset(x,x$month=="Sep-06")
Oct <- subset(x,x$month=="Oct-06")
Nov <- subset(x,x$month=="Nov-06")

    
 #define number of individuals to sample for each month
number <- c(80,77,84)  


#Setting the number of repetitions of the loop
n_reps <- 10000

#make 2 matrices to store the selection intensities for sept_oct and oct_nov
oct_nov <- sept_oct <- matrix(NA,nrow=n_reps,ncol=1)
colnames(oct_nov) <- colnames(sept_oct) <-c("intensity")

#Run the bootstrapped analysis
for(i in 1:n_reps){

		Sadjplates <- sample(x = Sept$adj.plates, size = 80, replace=TRUE, prob = NULL)
		Oadjplates <- sample(x = Oct$adj.plates, size = 77, replace=TRUE, prob = NULL)
		Nadjplates <- sample(x = Nov$adj.plates, size = 84, replace=TRUE, prob = NULL)
		
		
		means.sept <- mean(Sadjplates)
		sdS <- sd(Sadjplates)
		means.oct <- mean(Oadjplates)
		sdO <- sd(Oadjplates)
		means.nov <- mean(Nadjplates)
		sdN <- sd(Nadjplates)

		i_SO <- ((means.oct - means.sept)/sdS)
		sept_oct[i,] <- sept.oct <- unlist(i_SO)
		i_ON <- ((means.nov - means.oct)/sdO)
		oct_nov[i,] <- oct.nov <- unlist(i_ON)
}

#Outputting the results
mean(sept_oct[,1])
0.03929486
quantile(sept_oct[,1], probs=seq(0,1,0.025))
#range 2.5% to 97.5%
-0.292097571 - 0.381192429

mean(oct_nov[,1])
-0.03972177
quantile(oct_nov[,1], probs=seq(0,1,0.025))
-0.327927069 - 0.253828625


 
