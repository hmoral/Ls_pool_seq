#########################################
##    Cluster Separation Score (CSS)   ##
##    from an allele frequency file    ##
##    CSS meassure between-ecotype     ##
##    weighted by within-ecotype       ## 
##          genetic distances          ##
#########################################

# written by Hernán E. Morales
# if the paper is not published by the time of use please cite both the paper in press and the preprint
# paper: Morales, H. E., Faria, R., Johannesson, K., Larsson, T., Panova, M., Westram, A. M., & Butlin, R. (in press) Genomic architecture of parallel ecological divergence: beyond a single environmental contrast. Science Advances
# preprint: Morales, H. E., Faria, R., Johannesson, K., Larsson, T., Panova, M., Westram, A. M., & Butlin, R. (2018). Genomic architecture of parallel ecological divergence: beyond a single environmental contrast. bioRxiv, 447854.

# Changelog

# INPUT: 4 arguments
# The script takes two files (arguments 1 and 2)

# 1) The first, is the allelic frequencies file with SNPs in rows and pools in columns
# Allelic frequencies are from the reference allele
# Frequencies should laready be filtered by quality and minor allele frequency (downstream)
# The file must contain information about every SNP:
# Linkage Group or chromosome
# Position in cM within the LG or chrom
# The position on LG or chrom can be changed for contig or scaffold if desired to apply CSS to a smaller genomic scale (i.e. one CSS per contig/scaffold)
# example file for SNP frequencies for 6 pools

# LG  cM  Sw5_Crab   Sw5_Wave   Sw3_Crab   Sw3_Wave  Sw4_Crab  Sw4_Wave
# 1:  1  3.25566 0.3653846 0.15942029 0.44615385 0.22666667 0.2941176 0.4400000
# 2:  1  3.25566 0.1296296 0.11111111 0.20312500 0.17333333 0.1764706 0.2692308
# 3:  1  3.25566 0.3636364 0.20289855 0.44444444 0.26315789 0.2857143 0.4230769
# 4:  1  3.25566 0.1529412 0.13333333 0.21875000 0.15254237 0.2857143 0.1153846
# ---                                                                           
# 9581: 17 35.97442 0.1891892 0.16071429 0.23214286 0.06557377 0.1250000 0.3243243
# 9582: 17 35.97442 0.2027027 0.13207547 0.22641509 0.07575758 0.1111111 0.2631579
# 9583: 17 35.97442 0.2812500 0.22641509 0.24489796 0.44000000 0.2608696 0.1515152
# 9584: 17 35.97442 0.3492063 0.48837209 0.17333333 0.24719101 0.1176471 0.3478261

# 2) The second file, contains the information about the pools
# This is needed to differentiate genetic distances of between-ecotypes from those within-ecotypes
# example file for the frequency file above

# code ecotype
# 1: Sw3_Crab    Crab
# 2: Sw3_Wave    Wave
# 3: Sw4_Crab    Crab
# 4: Sw4_Wave    Wave
# 5: Sw5_Crab    Crab
# 6: Sw5_Wave    Wave

# 3) Number of pools included in the frequency file
# 4) Number of bootstrap repetitions
# 5) Prefix for output file names

# REQUIRED
# R-libraries data.table

# OUTPUT 
# A single table with observed and bootstrap CSS values (see at the bottom of script for example

args <- commandArgs(trailingOnly = TRUE)
## Default setting when no arguments passed
if(length(args) < 5) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      Arguments:
      --arg1= allelic frequencies file
      --arg2= sample infomation file
      --arg3= number of pools
      --arg4= number of bootstraps
      --arg5= output file prefix
      --help              - print this text
      
      Example:
      Rscript cluster_separation_score.R FREQ_FILE.txt.gz SampleInfo_FILE.txt.gz 200 out_CSS \n\n")
  
  q(save= "no")
}

##################
# Load libraries #
##################
library(data.table)

##################
# Load arguments #
##################
FILE<-args[[1]]
SampleInfo<-args[[2]]
poolsNo<-as.numeric(args[[3]])
bootstrap<-as.numeric(args[[4]])
OUTprefix<- args[[5]]

FILE<-"/Users/Fini/Documents/Littorina/Pool-Seq/1_FINAL/1_mns/1_scripts/test_trash.txt.gz"
SampleInfo<-"~/Documents/Littorina/Pool-Seq/1_FINAL/PCA/sampleInfo.txt"
poolsNo<-22
bootstrap<-10
OUTprefix<- "out_test"


##################
#   Script       #
##################

# read file in termination ".gz" is recognised as compressed file
if(grepl(".gz", FILE)){
  FREQfile<-fread(paste("zcat<",FILE),header=T)}else{FREQfile<-fread(FILE,header=T)
  }
# read file
SampleInfo<-fread(SampleInfo,h=T)

# generate a table where the bootstrap information will be stored
bootstrap_table<-data.table(matrix(data = NA,nrow = 1,ncol = bootstrap))
names(bootstrap_table)<-paste("b",1:bootstrap,sep="")
# generate a table where the observed values will be stored
# here we have two ecotypes (Crab and Wave) and columns represent genetic distances, their sd, and the CSS values
CSS_table<-data.table(LG=NA,cM=NA,SNPs=NA,eco_dist=NA,eco_dist_sd=NA,crab_dist=NA,crab_dist_sd=NA,wave_dist=NA,wave_dist_sd=NA,CSS=NA)
# combine tables
CSS_table<-cbind(CSS_table,bootstrap_table)
CSS_table<-CSS_table[0,]
COUNT<-0
PROG<-0
totalSNPS<-nrow(FREQfile)
# loop over LGs
for(i in unique(FREQfile$LG)){
  # print LG that is being processed and the progress in %
  print(paste("LG",i,"progress",PROG))
  # CSS scored are calculated for each unique cM per LG
  for(j in unique(FREQfile[LG==i]$cM)){
    COUNT<-COUNT+1
    PROG<-(COUNT*100)/totalSNPS
    FREQfile_section<-FREQfile[LG==i&cM==j]
    FREQfile_section<-FREQfile_section[,-c("LG","cM"),with=F]
    # PCA cannot be performed with less values than variables, here we have 22 pools, so position with less SNPs than 22 are discarded
    if(nrow(FREQfile_section)>poolsNo){
      mia.pca<-princomp(FREQfile_section)
      DIST_pca<-dist(mia.pca$loadings[,1:4],method = "euclidean")
      # generate a table with PCA distances for each ecotype 
      m <- data.table(t(combn(rownames(mia.pca$loadings),2)), as.numeric(DIST_pca))
      names(m) <- c("code", "pop2", "distance")
      m<-merge(m,SampleInfo[,c("code","ecotype"),with=F],by = "code")
      names(m)[grep("^code",names(m))]<-"pop1"
      names(m)[grep("ecotype",names(m))]<-"ecotype1"
      names(m)[grep("pop2",names(m))]<-"code"
      m<-merge(m,SampleInfo[,c("code","ecotype"),with=F],by = "code")
      names(m)[grep("^code",names(m))]<-"pop2"
      names(m)[grep("ecotype$",names(m))]<-"ecotype2"
      # distance calculations between and within ecotypes
      # formula from Supplementary Material of Jones et al., 2012; Nature, page 5 (https://doi.org/10.1038/nature10944).
      eco_distance=sum(m$distance[((m$ecotype1=="Crab"&m$ecotype2=="Wave")|
                                     (m$ecotype1=="Wave"&m$ecotype2=="Crab")|
                                     (m$ecotype2=="Crab"&m$ecotype1=="Wave")|
                                     (m$ecotype2=="Wave"&m$ecotype1=="Crab"))])
      crab_distance=sum(m$distance[(m$ecotype1=="Crab"&m$ecotype2=="Crab")])
      wave_distance=sum(m$distance[(m$pop1!=m$pop2&m$ecotype1=="Wave"&m$ecotype2=="Wave")])
      crab_sampleSize=poolsNo/2;wave_sampleSize=poolsNo/2
      ecotype_sum=eco_distance/(crab_sampleSize*wave_sampleSize)
      crab_sum=crab_distance/(crab_sampleSize^2*(crab_sampleSize-1))
      wave_sum=wave_distance/(wave_sampleSize^2*(wave_sampleSize-1))
      CSS = ecotype_sum-(crab_sampleSize+wave_sampleSize)*(crab_sum+wave_sum)
      CSS_all=CSS
      rm(list=c("CSS",ls(pattern = "_distance"),ls(pattern = "_sum")))
      distance_table_slice<-data.table(slice=NA,eco_dist=NA,crab_dist=NA,wave_dist=NA,CSS=NA)
      distance_table_slice<-distance_table_slice[0,]
      # repeat the same for X number of bootsrap
      # for each bootstrap the same number of SNPs are randomly sampled X number of times across the genome
      for(z in 1:bootstrap){
        FREQfile_section_slice<-FREQfile_section[sample(1:nrow(FREQfile_section),nrow(FREQfile_section),replace = T)]
        mia.pca<-princomp(FREQfile_section_slice)
        DIST_pca<-dist(mia.pca$loadings[,1:4],method = "euclidean")
        m <- data.table(t(combn(rownames(mia.pca$loadings),2)), as.numeric(DIST_pca))
        names(m) <- c("code", "pop2", "distance")
        m<-merge(m,SampleInfo[,c("code","ecotype"),with=F],by = "code")
        names(m)[grep("^code",names(m))]<-"pop1"
        names(m)[grep("ecotype",names(m))]<-"ecotype1"
        names(m)[grep("pop2",names(m))]<-"code"
        m<-merge(m,SampleInfo[,c("code","ecotype"),with=F],by = "code")
        names(m)[grep("^code",names(m))]<-"pop2"
        names(m)[grep("ecotype$",names(m))]<-"ecotype2"
        # distance calculations
        eco_distance=sum(m$distance[((m$ecotype1=="Crab"&m$ecotype2=="Wave")|
                                       (m$ecotype1=="Wave"&m$ecotype2=="Crab")|
                                       (m$ecotype2=="Crab"&m$ecotype1=="Wave")|
                                       (m$ecotype2=="Wave"&m$ecotype1=="Crab"))])
        crab_distance=sum(m$distance[(m$ecotype1=="Crab"&m$ecotype2=="Crab")])
        wave_distance=sum(m$distance[(m$pop1!=m$pop2&m$ecotype1=="Wave"&m$ecotype2=="Wave")])
        crab_sampleSize=11;wave_sampleSize=11
        # scores calculations
        ecotype_sum=eco_distance/(crab_sampleSize*wave_sampleSize)
        crab_sum=crab_distance/(crab_sampleSize^2*(crab_sampleSize-1))
        wave_sum=wave_distance/(wave_sampleSize^2*(wave_sampleSize-1))
        CSS = ecotype_sum-(crab_sampleSize+wave_sampleSize)*(crab_sum+wave_sum) # corrected by Roger
        distance_table_slice_x<-data.table(slice=z,eco_dist=eco_distance,crab_dist=crab_distance,wave_dist=wave_distance,CSS=CSS)
        distance_table_slice<-rbind(distance_table_slice,distance_table_slice_x)
      }
      # populate bootstrap_table for this step and combine with main table
      bootstrap_table_x<-data.table(matrix(data = distance_table_slice$CSS,nrow = 1,ncol = bootstrap))
      names(bootstrap_table_x)<-paste("b",1:bootstrap,sep="")
      # populate  CSS_table for this step and combine with main table
      CSS_table_x<-data.table(LG=i,cM=j,SNPs=nrow(FREQfile_section),eco_dist=mean(distance_table_slice$eco_dist),eco_dist_sd=sd(distance_table_slice$eco_dist),
                                       crab_dist=mean(distance_table_slice$crab_dist),crab_dist_sd=sd(distance_table_slice$crab_dist),
                                       wave_dist=mean(distance_table_slice$wave_dist),wave_dist_sd=sd(distance_table_slice$wave_dist),
                                       CSS=CSS_all)
      # combine both tables
      CSS_table_x<-cbind(CSS_table_x,bootstrap_table_x)
      CSS_table<-rbind(CSS_table,CSS_table_x)
    }
  }
}

# The output is a table that contains the following for each cM map positon: 
# mean genetic distances between and within ecotypes and their SD
# Number of SNPs contained in that cM map position
# CSS for observed values
# CSS for randomly selected loci b1,b2,b3... up to # of bootstraps
# example:
# LG        cM SNPs eco_dist eco_dist_sd crab_dist crab_dist_sd wave_dist wave_dist_sd        CSS         b1         b2         b3
# 1:  1  3.255660  166 49.54293   0.3054737  23.52306   0.38155295  24.21625   0.37813026 -0.4604713 -0.4545816 -0.4619812 -0.4612309
# 2:  1 26.178034  110 50.29149   0.7471649  22.48921   1.02222394  23.36382   0.64430686 -0.4172576 -0.4321428 -0.4019021 -0.4113936
# 3:  1 32.324249  106 51.20338   0.4274369  26.76310   1.02732953  20.87508   1.40959620 -0.4445326 -0.4427498 -0.4325469 -0.4395813
# 4:  1 35.088499  108 53.79648   0.4576694  25.64831   0.60608653  20.79122   0.57763337 -0.4079137 -0.4063317 -0.4014991 -0.3932820
# 5:  1 35.643113  177 54.94990   0.3028157  23.44105   0.35432366  22.40165   0.54197895 -0.3718685 -0.3801592 -0.3715977 -0.3808992

# a compressed table is output

gz1 <- gzfile(paste(OUTprefix,"_CSS_table_axes14_boot_euclidean.txt.gz",sep=""), "w")
write.table(CSS_table, gz1,row.names = F,quote = F,sep="\t")
close(gz1)