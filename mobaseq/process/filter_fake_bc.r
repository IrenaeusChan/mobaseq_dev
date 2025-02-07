#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

#setwd("/storage1/fs1/bolton/Active/Users/IrenaeusChan/MOBASeq-Toolkit/tests/MergeReads_2024_3DMobaLiv_ZY-CXCR2")
#setwd("/storage1/fs1/bolton/Active/Users/IrenaeusChan/MOBASeq-Toolkit/tests/")
setwd(args[1])

sink("proceedings_log.rout")

start_time <- Sys.time()
cat("Start time:", format(start_time), "\n")

#SampleFinal <- read.table(file = "Mobaseq_Sample_Info.txt", sep = "\t", stringsAsFactors = F, header = T)
SampleFinal <- read.table(file = args[2], sep = "\t", stringsAsFactors = F, header = T)

# 2Cluster tumor and remove spurious tumors.R
# @author: Chuan Li (chuanli@stanford.edu)

### Cluster reads and remove spurious tumors
# Any "tumours" that are within a Hamming distance of two from a larger tumour is assigned as "spurious tumours", 
# which are likely to be resulting from sequencing or PCR error, and are removed from subsequent analysis. 
# we only focus on the ones with the correct length 


## keep only the one that are more than 3 steps away from the focal sequence
myDist <- function(Str1){
  seq1<-unlist(strsplit(Str1,split=""))
  seq2<-unlist(strsplit(myBC,split=""))
  # check if they are the same length
  stopifnot(length(seq1) == length(seq2))
  Dist=0
  for (i in 1:length(seq1)){
    if (seq1[i] != seq2[i]){
      Dist = Dist + 1
    }
    if (Dist > 3){
      return(TRUE) # keep the sequence
    }
  }
  return(FALSE) # filter the sequence
}

myDistGB <- function(Str1){
  seq1<-unlist(strsplit(Str1,split=""))
  # check if they are the same length
  stopifnot(length(seq1) == length(myBC))
  return(as.logical( sum(seq1 != myBC) > 2  ))
}


for (i in 1:dim(SampleFinal)[1]){
  
  #FileNames <- paste0(SampleFinal$MergeReadOut[i], "_BarcodeClean.txt")
  # Data <- read.csv(file = SampleFinal$MergeReadOut[i], header = T, sep="\t", stringsAsFactors = F) # GB***: this is never used later, I commented it out
  
  #change filename and output file name
  #FileName <- paste0(SampleFinal$MergeReadOut[i], "_sgIDCounts.txt") # replace InputFile with the input
  #FileOut <- sub(pattern = "(*).csv_sgIDCounts.txt", "\\1_BarcodeClean.txt", FileName)

  FileName <- SampleFinal$MergeReadOut[i]
  FileOut <- sub("_MergeReadOut.csv$", "_BarcodeClean.txt", FileName)
  cat("\n\n\n----------- Starting ",FileName," -------------\n\n")
  
  #add "sgID","BC","counts"as column name
  SamplesFull <- read.csv(file = SampleFinal$MergeReadOut[i], header = F, sep=",", stringsAsFactors = F)
  names(SamplesFull) <- c("sgID","dist","BC","Count")
  SamplesFull <- SamplesFull[order(-SamplesFull$Count),]
  SamplesFull <- SamplesFull[SamplesFull$Count>1,] # GB: onle BCs with at least 2 reads are kept
  
  #empty vector
  SamplesClean <- c()
  
  #all the sgID in MergeReadOut
  sgIDAll <- unique(SamplesFull$sgID)
  
  cat("sgIDs to loop through: ",sgIDAll)
  # find whether the sequence contain any N
  # remove these ones
  findN <- function(String){
    return(grepl('N',String))
  }
  
  BoolfindN <- sapply(SamplesFull$BC, FUN = findN)
  SamplesFull <- SamplesFull[!BoolfindN, ]
  
  
  for (sgID in sgIDAll){ 
    cat("\nStarting to filter: ",sgID,"\n")
    Samples <- SamplesFull[SamplesFull$sgID == sgID,]
    N <- dim(Samples)[1]
    BoolAll <- rep(TRUE, N)
    
    
    cat("Number of barcodes to test: ", N,"\n\n")
    for (iBC in 1:(N-1)){  # GB: N-1 because for the last one, there is nothing else to compare it to.
      if (iBC %% 1000 == 0){ cat(iBC,"_");}
      
      if (Samples$Count[iBC] < 5){ # GB: if a BC has 4 or less reads, we are not comparing it to any smaller BCs
        break
      }
      
      if(N == 1) { # AX240131: if a BC has no other valid BCs (reads > 1) to compare to, move on to the next one to prevent a bug 
        break
      }
      
      # For current barcode
      # myBC <- Samples$BC[iBC]
      myBC <- unlist(strsplit(Samples$BC[iBC],""))
      
      # Filter through the rest of barcode
      
      
      if (iBC == N-1){
        DistAll <- myDistGB(Samples$BC[(iBC+1):N])
      }else{
        DistAll <- sapply(Samples$BC[(iBC+1):N], FUN = myDistGB)
      }
      
      DistAll <- c(rep(TRUE, iBC), DistAll)
      BoolAll[!DistAll] <- FALSE
    }
    
    Samples <- Samples[BoolAll, ]
    cat("\nPercent of BCs filtered:", signif(1 - mean(BoolAll) , 2)*100, "%\n\n")
    SamplesClean <- rbind(SamplesClean, Samples)
    # GB: write a temporary file
    write.table(SamplesClean, file = sub(".txt","_temporary",FileOut), append = FALSE, quote = FALSE, sep = "\t", row.names = F,
                col.names = TRUE)
  }
  
  write.table(SamplesClean, file = FileOut, append = FALSE, quote = FALSE, sep = "\t", row.names = F,
              col.names = TRUE)
  file.remove(sub(".txt","_temporary", FileOut))
}

end_time <- Sys.time()
cat("End time:", format(end_time), "\n")

duration <- end_time - start_time
cat("Duration:", round(duration, 2), "seconds\n")

sink()