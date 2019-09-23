#########################################
## Classification of depth into three  ##
## classes of coverage distributions:  ##
##    low, medium, and high            ##
## with a mixed distribution model     ##
#########################################

# written by Hernán E. Morales
# if the paper is not published by the time of use please cite both the paper in press and the preprint
# paper: Morales, H. E., Faria, R., Johannesson, K., Larsson, T., Panova, M., Westram, A. M., & Butlin, R. (in press) Genomic architecture of parallel ecological divergence: beyond a single environmental contrast. Science Advances
# preprint: Morales, H. E., Faria, R., Johannesson, K., Larsson, T., Panova, M., Westram, A. M., & Butlin, R. (2018). Genomic architecture of parallel ecological divergence: beyond a single environmental contrast. bioRxiv, 447854.

# Changelog

# INPUT: 4 arguments
# 1) The script takes one file with position in the row (500 bp in the example below) and at least one column with a coverage value (cov_mean in the example below)
# Note: other formats for the depth metrics should work

#       contig chromStart chromEnd     cov_mean      cov_sd
# 1:      Contig0          0      500 15.490000000  7.53637515
# 2:      Contig0        500     1000 70.767333333 33.69770994
# 3:      Contig0       1000     1500  0.471000000  0.49444231
# 4:      Contig0       1500     2000 60.693666667 11.52212851
# 5:      Contig0       2000     2500  0.245000000  0.17719029
# ---                                                          
# 3202654: Contig388609       1500     2000 18.999333333  4.30394157
# 3202655: Contig388612          0      500  0.006666667  0.01632993
# 3202656: Contig388612        500     1000  8.459000000  1.76396542
# 3202657: Contig388613          0      500  0.000000000  0.00000000
# 3202658: Contig388613        500     1000 25.692000000  3.20829849

# 2) The name of the column that contains the coverage values
# 3) Number of iterations of random subsampling
# 4) prefix for output file names


# REQUIRED
# R-libraries data.table, mixtools

# OUTPUT 
# a plot of the raw distribution and the mixed model outputs for each iteration
# a single table with the limits of the 50% distribution for each coverage class

args <- commandArgs(trailingOnly = TRUE)
## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      Arguments:
      --arg1= coverage depth file
      --arg2= column name for coverage
      --arg3= number of iterations
      --arg4= output file prefix
      --help              - print this text
      
      Example:
      Rscript mixed_dist_model_depth.R FILE.txt.gz cov_mean 200 out_mixMod \n\n")
  
  q(save= "no")
}

##################
# Load libraries #
##################
library(data.table)
library(mixtools)

##################
# Load arguments #
##################
FILE<- args[[1]]
COV<- args[[2]]
ITER<- as.numeric(args[[3]])
OUTprefix<- args[[4]]

##################
#   Script       #
##################
# read file in termination ".gz" is recognised as compressed file
if(grepl(".gz", FILE)){
  depth_file<-fread(paste("zcat<",FILE),header=T)}else{depth_file<-fread(FILE,header=T)
}
# define coverage column
COVcolumn<-grep(COV,names(depth_file))
names(depth_file)[COVcolumn]<-"COVERAGE"

# store coverage values in a vector, remove cero values and NAs
DP_COVERAGE<-depth_file$COVERAGE[depth_file$COVERAGE>0&!is.na(depth_file$COVERAGE)]

# The mixmdl becaomes slow when more than 500K records are analyzed. 
# One solution is to iterate over the data and randomly sample 500K values in each iteration. 
# This is implemented here with a simple for loop

# Create a table that will receive the value in each iteration
TABLE<-data.frame(Iter=NA,K1.L1=NA,K1.L2=NA,K2.L1=NA,K2.L2=NA,K3.L1=NA,K3.L2=NA)
TABLE<-TABLE[0,]
# start for loop
for(i in 1:ITER){
  OUT<-paste(OUTprefix,"_mixDist_iter_",i,sep="")
  DP_COVERAGE_sample<-sample(x = DP_COVERAGE,size = 500000,replace = F)
  png(paste(OUT,"_plots.png",sep=""))
  par(mfrow=c(3,2))
  plot(density(DP_COVERAGE_sample),xlim=c(0,300))
  DP.sqrt<-sqrt(DP_COVERAGE_sample)
  #DP.sqrt<-sample(DP.sqrt,50000,replace = F)
  mixmdl = normalmixEM(DP.sqrt,k = 3,epsilon=1e-3)
  plot(mixmdl,which=2)
  lines(density(DP.sqrt), lty=2, lwd=2)
  #summary(mixmdl)
  MEANS<-mixmdl$mu^2
  SDs<-mixmdl$sigma^2
  LMBDA<-mixmdl$lambda^2
  plot(0,xlim=range(mixmdl$mu^2),ylim=range(mixmdl$sigma^2),type="n",
       xlab="Component mean", ylab="Component standard deviation")
  points(x=mixmdl$mu^2,y=mixmdl$sigma^2,pch=as.character(1:6),
         cex=sqrt(0.5+5*(mixmdl$lambda^2)))
  #summary(mixmdl) # this is the result of the mixture model
  lambda1 <- mixmdl$lambda[1]
  lambda2 <- mixmdl$lambda[2]
  lambda3 <- mixmdl$lambda[3]
  mu1 <- mixmdl$mu[1]
  mu2 <- mixmdl$mu[2]
  mu3 <- mixmdl$mu[3]
  sigma1 <- mixmdl$sigma[1]
  sigma2 <- mixmdl$sigma[2]
  sigma3 <- mixmdl$sigma[3]
  
  x <- sort(round(DP_COVERAGE_sample,0))
  x<-unique(x)
  plot(density(x),main="density raw unique values")
  sqrt_x <- x^0.5
  plot(density(sqrt_x),main="density sqrt unique values")
  
  d1 <- (pnorm(sqrt_x[-1],mu1,sigma1)-pnorm(sqrt_x[-length(x)],mu1,sigma1))*lambda1
  #plot(x[-1],d1)
  d2 <- (pnorm(sqrt_x[-1],mu2,sigma2)-pnorm(sqrt_x[-length(x)],mu2,sigma2))*lambda2
  #points(x[-1],d2,col="blue")
  d3 <- (pnorm(sqrt_x[-1],mu3,sigma3)-pnorm(sqrt_x[-length(x)],mu3,sigma3))*lambda3
  #points(x[-1],d3,col="red")
  
  P_d1 <- d1/(d1+d2+d3)
  P_d2 <- d2/(d1+d2+d3)
  P_d3 <- d3/(d1+d2+d3)
  
  y <- x[-1]
  CUT_OFF<-0.5
  K1.1<-min(y[P_d1>CUT_OFF],na.rm = T)
  K1.2<-max(y[P_d1>CUT_OFF],na.rm = T)
  K2.1<-min(y[P_d2>CUT_OFF],na.rm = T)
  K2.2<-max(y[P_d2>CUT_OFF],na.rm = T)
  K3.1<-min(y[P_d3>CUT_OFF],na.rm = T)
  K3.2<-max(y[P_d3>CUT_OFF],na.rm = T)
  LIMTS<-sort(c(K1.1,K1.2,K2.1,K2.2,K3.1,K3.2))
  TABLEx<-data.frame(Iter=i,K1.L1=LIMTS[1],K1.L2=LIMTS[2],K2.L1=LIMTS[3],K2.L2=LIMTS[4],K3.L1=LIMTS[5],K3.L2=LIMTS[6])
  TABLE<-data.frame(rbind(TABLE,TABLEx))
  TITLE<-paste("K1=",LIMTS[1],LIMTS[2],"- K2=",LIMTS[3],LIMTS[4],"- K3=",LIMTS[5],LIMTS[6])
  plot(x[-1],P_d1,xlim=c(0,500),main=TITLE)
  points(x[-1],P_d2,col="blue")
  points(x[-1],P_d3,col="red")
  # 
  dev.off()
}
write.table(TABLE.1,paste(OUT,"_iterations_table.png",sep=""),row.names = F,quote = F,sep="\t")
