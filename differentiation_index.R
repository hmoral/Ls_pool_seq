#########################################
##      Differentiation Index (DI)     ##
##    meassure the directionality of   ##
##  allelic frequency differentiation  ##
#########################################

# written by Hernán E. Morales
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
# example file for Low vs High pools (see Morales et al)

# LG     cM    Fr_Low    Fr_High  SPs_High    SPs_Low  SPn_High   SPn_Low SFreqTablen5_C_High SFreqTablen5_C_Low SFreqTablen5_FreqTable_High SFreqTablen5_FreqTable_Low SFreqTablen3_C_High SFreqTablen3_C_Low SFreqTablen3_FreqTable_High SFreqTablen3_FreqTable_Low  Ukw_High   Ukw_Low  Uke_High   Uke_Low
# 1:  1  0.657 0.6363636 0.75675676 0.4255319 0.75000000 0.7529412 0.6511628   0.7708333  0.5348837   0.2758621  0.3333333  0.80882353 0.63461538   0.5652174  0.3896104 0.5058824 0.6708861 0.5588235 0.5660377
# 2:  1  0.657 0.6666667 0.75652174 0.3958333 0.73076923 0.7790698 0.6595745   0.7307692  0.5581395   0.2941176  0.3333333  0.82608696 0.67272727   0.3809524  0.3780488 0.5308642 0.6585366 0.5514019 0.6160714
# 3:  1  0.657 0.1891892 0.06818182 0.1102362 0.14754098 0.1844660 0.3095238   0.1774194  0.3296703   0.2500000  0.2321429  0.05882353 0.07843137   0.1923077  0.1666667 0.3982301 0.2087912 0.1487603 0.1848739
# 4:  1  0.657 0.5304348 0.49618321 0.5632184 0.92857143 0.4122807 0.5000000   0.3278689  0.4337349   0.7959184  0.6833333  0.31250000 0.38095238   0.5000000  0.6428571 0.8888889 0.7727273 0.6564885 0.5816327
# ---                                                                                                                                                                                                            
# 9997:  1 12.850 0.3504274 0.39230769 0.8349515 0.73737374 0.5754717 0.8974359   0.8679245  0.8571429   0.7755102  0.7021277  0.92000000 0.93750000   0.9565217  0.7717391 0.8469388 0.7941176 0.7591241 0.7090909
# 9998:  1 12.850 0.4000000 0.53488372 0.8800000 0.74038462 0.5128205 0.8611111   0.9107143  0.8873239   0.7826087  0.7272727  0.93750000 0.82051282   0.9615385  0.8085106 0.7289720 0.7826087 0.8181818 0.8073394
# 9999:  1 12.850 0.3238095 0.20161290 0.4421053 0.15686275 0.1698113 0.2682927   0.7857143  0.7164179   0.6585366  0.4705882  0.42105263 0.55813953   0.2500000  0.4421053 0.3203883 0.4179104 0.4587156 0.4392523

# 2) The second file, contains the information about the pools
# This is needed to differentiate the differentiation index in each pool and to assign the reference pool to confer directionality
# example file for the frequency file above
# code shore country
# 1:      Fr_Low   Low  France
# 2:     Fr_High  High  France
# 3:    SPn_High  High   Spain
# 4:     SPn_Low   Low   Spain
# ---                                                                           
# 16:    Uke_High  High      UK
# 17:     Ukw_Low   Low      UK
# 18:    Ukw_High  High      UK

# 3) Reference group (country in our case): the group of reference used as reference to confer directionality to the Differentiation Index. Here this is assigned by defining the name of country
# 4) Prefix for output file names

# REQUIRED
# R-libraries data.table

# OUTPUT 
# A single table with DI values (see at the bottom of script for example)

args <- commandArgs(trailingOnly = TRUE)
## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      Arguments:
      --arg1= allelic frequencies file
      --arg2= sample infomation file
      --arg3= reference group
      --arg4= output file prefix
      --help              - print this text
      
      Example:
      Rscript differentiation_index.R FREQ_FILE.txt.gz SampleInfo_FILE.txt.gz Spain out_DI \n\n")
  
  q(save= "no")
}

##################
# Load arguments #
##################
FILE<-args[[1]]
SampleInfo<-args[[2]]
REFcountry<-args[[3]]
OUTprefix<- args[[4]]

##################
# Load libraries #
##################
library(data.table)

##################
#   Script       #
##################

# read frequency file in termination ".gz" is recognised as compressed file
if(grepl(".gz", FILE)){
  FreqTable<-fread(paste("zcat<",FILE),header=T)}else{FreqTable<-fread(FILE,header=T)
  }

# read file
SampleInfo<-fread(SampleInfo,h=T)

# define sample sets
ALL_pools<-as.character(SampleInfo$code)
HIGH_pools<-sort(as.character(SampleInfo[shore=="High"]$code))
LOFreqTable_pools<-sort(ALL_pools[!ALL_pools%in%HIGH_pools])

# define reference sets
high_REF<-as.character(SampleInfo[country==REFcountry&shore=="High"]$code)
low_REF<-as.character(SampleInfo[country==REFcountry&shore=="Low"]$code)

# make a new matrix to store results
DI_table<-data.frame(LG=NA,cM=NA,REF_mean=NA,matrix(data = NA,nrow = 1,ncol = length(HIGH_pools)))
names(DI_table)[4:ncol(DI_table)]<-HIGH_pools
DI_table<-DI_table[0,]
# Hybrid Index = HI
i=1
LG<-unique(FreqTable$LG)
# loop over LGs
for(i in LG){
  # susbset to LG
  FreqTable_LG<-FreqTable[LG==i]
  # loop over map positions
  cM_position<-unique(FreqTable_LG$cM)
  for(j in cM_position){
    # susbset to a map position
    FreqTable_LG_cM<-FreqTable_LG[cM==j]
    # take average for high shore level reference pools 
    HIGH_pools_table_average<-rowMeans(FreqTable_LG_cM[,high_REF,with=F])
    # take average for low shore level reference pools 
    LOFreqTable_pools_table_average<-rowMeans(FreqTable_LG_cM[,low_REF,with=F])
    # get allelic frequency differential
    DIFF<-HIGH_pools_table_average-LOFreqTable_pools_table_average
    # ask which positions have negative differences 
    POS<-which(DIFF<0)
    # then for those positions take the frequency of the other pool (i.e. 1-p)
    if(length(POS)>0){HIGH_pools_table_average[POS]<-1-HIGH_pools_table_average[POS]}
    if(length(POS)>0){LOFreqTable_pools_table_average[POS]<-1-LOFreqTable_pools_table_average[POS]}
    HIGH_pools_table_average<-rowMeans(FreqTable_LG_cM[,HIGH_pools,with=F])
    LOFreqTable_pools_table_average<-rowMeans(FreqTable_LG_cM[,LOFreqTable_pools,with=F])
    if(length(POS)>0){HIGH_pools_table_average[POS]<-1-HIGH_pools_table_average[POS]}
    if(length(POS)>0){LOFreqTable_pools_table_average[POS]<-1-LOFreqTable_pools_table_average[POS]}
    DI_table_x<-data.frame(LG=i,cM=j,REF_mean=(mean(HIGH_pools_table_average)-mean(LOFreqTable_pools_table_average)),
                                matrix(data = NA,nrow = 1,ncol = length(HIGH_pools)))
    names(DI_table_x)[4:ncol(DI_table_x)]<-HIGH_pools
    # then for each of the high pools, do the flipping and add its average DI to the table
    HIGH_pools_table<-FreqTable_LG_cM[,HIGH_pools,with=F]
    LOFreqTable_pools_table<-FreqTable_LG_cM[,LOFreqTable_pools,with=F]
    for(z in 1:length(HIGH_pools)){
      HIGH_pool_s<-HIGH_pools[z]
      LOFreqTable_pool_s<-LOFreqTable_pools[z]
      REPLACE<-as.numeric(unlist(1-HIGH_pools_table[POS,z,with=FALSE]))
      HIGH_pools_table[POS,names(HIGH_pools_table)[z] := REPLACE]
      REPLACE<-as.numeric(unlist(1-LOFreqTable_pools_table[POS,z,with=FALSE]))
      LOFreqTable_pools_table[POS,names(LOFreqTable_pools_table)[z] := REPLACE]
      COL2<-which(names(DI_table_x)==HIGH_pool_s)
      DI_table_x[,COL2]<- (mean(unlist(HIGH_pools_table[,z,with=F]))-mean(unlist(LOFreqTable_pools_table[,z,with=F])))
    }
    # merge each newly generated line with main table
    DI_table<-rbind(DI_table,DI_table_x)
  }
}

gz1 <- gzfile(paste(OUTprefix,"_DI_table_REF_",REFcountry,".txt.gz",sep=""), "w")
write.table(DI_table, gz1,row.names = F,quote = F,sep="\t")
close(gz1)

# example of output, the names of the pools are labelled "High" but they should be simply the name of the pool. This extra labell is a left-over from their original names in the Frequency Table
#    LG      cM   REF_mean       Fr_High   SPn_High   SPs_High   SWn3_C_High   SWn3_W_High   SWn5_C_High  SWn5_W_High     Uke_High      Ukw_High
# 1:  1  0.6570 0.04633239  0.0592410839 0.06822184 0.19828067 -3.116543e-04  3.608055e-02  0.0238851079  0.003072025  0.042658124 -0.0141362554
# 2:  1  0.9655 0.02852787  0.0210745941 0.12068976 0.08160229  2.194858e-02  3.484261e-05  0.0216379972  0.000129408 -0.006754697 -0.0036119589
# 3:  1  2.0905 0.02479449  0.0076523380 0.13540482 0.08602227  1.828682e-02 -1.036285e-02  0.0050609639 -0.006713134 -0.013195334  0.0009945203
# 4:  1  2.6565 0.02343077 -0.0078660603 0.10943387 0.06815485  1.894425e-02  6.462131e-03  0.0004925125  0.006616356 -0.007488625  0.0161276361
# ----
