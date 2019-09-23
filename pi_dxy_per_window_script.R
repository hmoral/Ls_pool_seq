############################################
#                                          #
# Script to calculate different components #
# of nucleotide diversity (pi) and Dxy     #
# from sync files of pool-seq data         #
#                                          #
############################################

# based on Gauthier et al. 2018; Mol. Ecol
# https://doi.org/10.5061/dryad.r9q3b
# adapted by Hernan E. Morales
# if the paper is not published by the time of use please cite both the paper in press and the preprint
# paper: Morales, H. E., Faria, R., Johannesson, K., Larsson, T., Panova, M., Westram, A. M., & Butlin, R. (in press) Genomic architecture of parallel ecological divergence: beyond a single environmental contrast. Science Advances
# preprint: Morales, H. E., Faria, R., Johannesson, K., Larsson, T., Panova, M., Westram, A. M., & Butlin, R. (2018). Genomic architecture of parallel ecological divergence: beyond a single environmental contrast. bioRxiv, 447854.

args <- commandArgs(trailingOnly = TRUE)
## Default setting when no arguments passed
if(length(args) < 6) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      R script to calcuate pi and dxy from sync files
      
      Arguments:
      --arg1= sync file input (ending with .sync or .sync.gz if compressed)
      --arg2= population A number in sync file where ncol(sync)-3 = pop1
      --arg3= population B number in sync file where ncol(sync)-3 = pop1
      --arg4= minimum depth
      --arg5= maximum depth
      --arg6= window size in bp
      --arg7= minimum read count (default = 2)
      --help              - print this text
      
      Example:
      Rscript pi_dxy_per_window_script.R FILE.sync.gz 7 8 20 240 1000 2
      
      Output:
      A table with the following columns:
      contig = name of contig or chromosomes
      start = window start coordinates
      end = window end coordinates
      length = length of genomic window
      sites = number of sites considered in the window
      sites_info = number of sites that produced output
      filt_percent = proportion of sites with NA values (i.e. filtered out) 
      PiA = pi in population A
      PiB = pi in population B
      PiT = pi total
      PiWth = average pi within popA and popB
      PiBtw = pi between popA and popB 
      Dxy = Dxy between popA and popB
      \n\n")
  
  q(save= "no")
}

# Load libraries
library(data.table)
library(stringr)

# Load arguments
syncFile_input<- args[[1]]
popA<- as.numeric(args[[2]])
popB<- as.numeric(args[[3]])
minDepth<-as.numeric(args[[4]])
maxDepth<-as.numeric(args[[5]])
WinSize<-as.numeric(args[[6]])
minCount<-as.numeric(args[[7]])
if(!exists("minCount")){minCount<-2}

# Generate output name
# the default output name is: XXX_pi_dxy_popXvspopX.txt.gz
if(str_sub(syncFile_input,-3,-1)==".gz"){
  Out_name<-gsub(".sync.gz","",syncFile_input)
} else {
  Out_name<-gsub(".sync","",syncFile_input)
}
Out_name<-str_split_fixed(Out_name,"/",Inf);Out_name<-Out_name[length(Out_name)]
Out_name<-paste(Out_name,"_window_pi_dxy_","pop",popA,"vs","pop",popB,".txt.gz",sep="")



# pi_dxy_function
# this function loops over rows in the genomic window
# calculates pi per pool, pi total, pi within, pi between and Dxy
# average metrics per window
# outputs per window:
# number of sites
# number of sites that generated output
# proportion of sites with NA values (i.e. filtered out) 
# metrics pi per pool, pi within, pi between and Dxy
# invariant sites are taken into account

Pi_out_list = c()
#pop1<-7;pop2<-8
pi_dxy_function <- function(pop1,pop2){
  assign(paste("table", pop1,pop2, sep=""), c())
  # loop over rows
  for (i in 1:slice){
    assign(paste("row", i, sep = ""), rbind(as.numeric(strsplit(as.character(pop_data[i,pop1]), ":")[[1]])[1:4],
                                             as.numeric(strsplit(as.character(pop_data[i,pop2]), ":")[[1]])[1:4]))
    # collect counts from sync file
    countsA= c(get(paste("row", i, sep = ""))[1,1] , get(paste("row", i, sep = ""))[1,2] , get(paste("row", i, sep = ""))[1,3] , get(paste("row", i, sep = ""))[1,4])
    countsB= c(get(paste("row", i, sep = ""))[2,1] , get(paste("row", i, sep = ""))[2,2] , get(paste("row", i, sep = ""))[2,3] , get(paste("row", i, sep = ""))[2,4])
    # if count == NA asign a huge value that will be eventually filtered out
    countsA[is.na(countsA)]<-1000
    countsB[is.na(countsB)]<-1000
    # count how many alleles there are, if Alleles > 2 they will be filtered out later on
    AllelesA<-sum(countsA>0);AllelesB<-sum(countsB>0)
    # collect metrics for pi estimates
    nA = sum(countsA)
    nB = sum(countsB)
    countsT = countsA + countsB
    nT = sum(countsT)
    Nuc = which.max(countsA)
    freqA = countsA[Nuc]/nA
    freqB = countsB[Nuc]/nB
    freqT = (countsA[Nuc] + countsB[Nuc])/nT
    # estimate pi for popA, popB and pi total
    PiA = (nA/(nA-1))*(1 - ((countsA[1])/nA)^2 - ((countsA[2])/nA)^2 - ((countsA[3])/nA)^2 - ((countsA[4])/nA)^2)
    PiB = (nB/(nB-1))*(1 - ((countsB[1])/nB)^2 - ((countsB[2])/nB)^2 - ((countsB[3])/nB)^2 - ((countsB[4])/nB)^2)
    # pi total is calculated with frequencies to avoid bias in positions where the read counts are heavily unbalanced between pools 
    PiT = (nT/(nT-1))*(1 - (((countsA[1]/nA)+(countsB[1]/nB))/2)^2 - (((countsA[2]/nA)+(countsB[2]/nB))/2)^2 - (((countsA[3]/nA)+(countsB[3]/nB))/2)^2 - (((countsA[4]/nA)+(countsB[4]/nB))/2)^2)
    # if pi values were estimated for both pools then calculate pi_within, pi_between and Dxy
    # pi per pool from Tajima 1989; Genetics. 123 (3): 585–95.
    # pi within and between adapted from Charlesworth et al. 1997; Genet. Res. 70, 155–174.
    # Dxy formula from Nei 1987 Columbia university press, eq 10.20
    if(!is.na(PiA)&!is.na(PiB)){PiWth = (PiA+PiB)/2; PiBtw = PiT - PiWth}
    if(!is.na(freqA)&!is.na(freqB)){Dxy = (freqA * (1-freqB)) + (freqB * (1-freqA))}
    # if filters are not meet asign NA to all output metrics
    if((nA>0&nA<=minCount)|((nB>0&nB<=minCount))){PiA = NA;PiB = NA;PiT = NA;PiWth = NA;PiBtw = NA;Dxy = NA}
    if((nA<minDepth|nA>maxDepth)|((nB<minDepth|nB>maxDepth))){PiA = NA;PiB = NA;PiT = NA;PiWth = NA;PiBtw = NA;Dxy = NA}
    if((AllelesA>2)|(AllelesB>2)){PiA = NA;PiB = NA;PiT = NA;PiWth = NA;PiBtw = NA;Dxy = NA}
    # update list with output metrics
    assign(paste("table", pop1,pop2, sep = ""), rbind(get(paste("table", pop1,pop2, sep = "")), c(PiA,PiB,PiT,PiWth,PiBtw,Dxy)))
  }
  # collect pi metrics and average over window
  pi_metrics = get(paste("table", pop1,pop2, sep=""))
  sites= length(pi_metrics[,1])
  NAcount = sum(is.na(pi_metrics[,1]))
  m_Pi = sum(!is.na(pi_metrics[,1]))
  if(m_Pi>0){NAcount = NAcount/sites}
  if(m_Pi==0){NAcount = 1}
  if(sum(is.na(pi_metrics[,1]))==length(pi_metrics[,1])){
    Pi_out_list = NA
  } else{
    assign(paste("window", pop1, sep=""), colSums(get(paste("table", pop1,pop2, sep="")), na.rm = T)/m_Pi)
    Pi_out_list = get(paste("window", pop1, sep=""))
  }
  # collect final output
  final_pi_metrics = c(sites,m_Pi,NAcount,Pi_out_list)
  return(final_pi_metrics)
}


# read in sync file
print(syncFile_input)
#syncFile = fread(cmd = paste("zcat",syncFile_input,"| head -100000"), sep = "\t")
syncFile = fread(syncFile_input, sep = "\t")
# add header
NoPops<-ncol(syncFile)-3
#syncFile[1:5,1:5]
names(syncFile)<-c("contig","position","ref",paste("pop",1:NoPops,sep=""))
# print info
print(paste("No. of pops in sync file = ",NoPops))
print(paste("pops selected = ",popA,"--",popB))
# generate gen list
contigList<-unique(syncFile$contig)
# generate progress bar
z<-0
imax<-length(contigList)
pb <- txtProgressBar(min = 0, max = imax, style = 3)
# loop over sync file, generate windows and calculate metrics
Pi_results<-data.table()
for(j in contigList){
  tryCatch({
    z<-z+1
    syncFile_contig<-syncFile[contig==j]
    syncFile_contig_lenght<-seq(1:nrow(syncFile_contig))
    SLICES<-split(syncFile_contig_lenght, ceiling(seq_along(syncFile_contig_lenght)/WinSize))
    for(i in SLICES){
      START<-i[1]
      END<-i[length(i)]
      syncFile_win = syncFile_contig[START:END]
      pop_data = data.frame(syncFile_win[,4:ncol(syncFile_win)])
      slice = nrow(pop_data)
      launch_function= pi_dxy_function(popA,popB)
      Pi_results_x<-data.table(contig=j,start=START,end=END,length=(END-START)+1,
                               sites=launch_function[1],sites_info=launch_function[2],filt_percent=launch_function[3],
                               PiA=launch_function[4],PiB=launch_function[5],PiT=launch_function[6],PiWth = launch_function[7],PiBtw=launch_function[8],Dxy=launch_function[9])
      Pi_results<-rbind(Pi_results,Pi_results_x)
    }
    setTxtProgressBar(pb, z)
  }, error=function(e){})
}
gz1 <- gzfile(Out_name, "w")
write.table(Pi_results, gz1,row.names = F,quote = F,sep = "\t")
close(gz1)
