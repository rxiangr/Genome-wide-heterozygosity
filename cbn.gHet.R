#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("3 argument must be supplied, \nUsage: Rscript --vanilla cbn.gHet.R <../yourpathOf*fe.txt>  <../yourpath of output> <your outputname>", call.=FALSE)
}

library(data.table)

#---read in some arguments
chrfepath <- args[1] #(e.g, hetfeout)
outputpath <- args[2] #(e.g., '.')
outputfn <- args[3] #(e.g., 'test')
dir.create(outputpath)

#---create file list
fefilelist <- list.files(chrfepath,'fe.txt')
snpnlist <- list.files(chrfepath,'fe.nSNPs')
#---read in chromosome-wise fe
felist <- list()
for (i in fefilelist)
{nobs <- nrow(fread(paste0(chrfepath,'/',fefilelist[1])))
 if(nrow(nobs==fread(paste0(chrfepath,'/',i))))
 {felist[[i]] <- fread(paste0(chrfepath,'/',i))[,-1]
 }else {stop('Number of rows (individuals) differ in some chromosome fe files') }
}
fedt <- do.call(cbind,felist)
csumdt <- data.frame(fread(paste0(chrfepath,'/',i))[,-2],rowSums(fedt))
#---read on snp number
loopsnpnlist <- list()
for (j in snpnlist)
 {loopsnpnlist[[j]] <- fread(paste0(chrfepath,'/',j))
 }
nSNPs <- sum(unlist(loopsnpnlist))
#---combine results
csumdt[,2] <- cbind(csumdt[,1],genomean=csumdt[,2]/nSNPs)
write.table(csumdt,paste0(outputpath,'/',outputfn,'.cbfe.qc.txt'),row.names=F,col.names=F,quote=F,sep='\t')

cat(paste0('Combining FE finished at ',Sys.time()),sep='\n')
cat(paste0('Results saved to ',outputpath,'/',outputfn,'.cbfe.qc.txt'),sep='\n')
