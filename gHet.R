#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("4 argument must be supplied, \nUsage: Rscript gHet.R <../yourpath/prefix of plink file> <../yourpath/full name of plink file.frq> <../yourpath of output> <your outputname>", call.=FALSE)
}

library(BEDMatrix)
library(data.table)

#---read in some arguments
genphathpref <- args[1] #(e.g, testplink/test.chr24)
afpathfn <- args[2] #(e.g., testplink/test.chr24.frq)
outputpath <- args[3] #(e.g., 'hetfeout')
outpref <- args[4] #(e.g., 'chr24')
dir.create(outputpath)

#--read in plink binary
gmt <- BEDMatrix(paste0(genphathpref,'.bed'))
#--read in bed file
afdt.ref <- fread(afpathfn)
#---check SNP order matach
if (!sum(sub("\\_.*", "",colnames(gmt))==unlist(afdt.ref$SNP))==ncol(gmt)) {stop('Please supply AF file in which the SNP name match with the *.bim file')}
#--make sure af is calculated on reference allele, A1 in this case
afdt.ref[,ref.af:=ifelse(A1==1,MAF,ifelse(A1==2,1-MAF,NA))]
af.vec <- as.numeric(unlist(afdt.ref[,7]))

#---recode genotype base on af
code.homoref <- round(-2*((1-af.vec)^2),3)
code.het <- round(2*(1-af.vec)*af.vec,3)
code.homoalt <- round(-2*(af.vec^2),3)
#---updated genotype matrix using block (10k size) approach
blocksize <- 5000
blckvec <- sort(unique(c(seq(0,ncol(gmt), by = blocksize),ncol(gmt))))
totNblocks <- length(blckvec)-1
cat(paste0('Loop for block-wise analysis started at ',Sys.time()),sep='\n')
cat(paste0('Total number of SNPs: ',ncol(gmt)),sep='\n')
cat(paste0('Size of each block: ',blocksize),sep='\n')
cat(paste0('Total number of blocks: ',totNblocks),sep='\n')
#--loop for blocks(memory efficient)
blocklist <- list()
#for (blckn in (1:5))
for (blckn in (1:totNblocks))
{blckendst <- blckvec[blckn]+1
 blckendend <- blckvec[blckn+1]
 test <- setDT(data.frame(gmt[,blckendst:blckendend]))
 test[] <- lapply(test, function(x) as.numeric(x))
 snplist <- colnames(test)
 reportlist <- seq(1,length(snplist),1000)
 #j <- 1
  for (j in 1:length(snplist))
   {if(j %in% reportlist) {cat(paste0('analysing SNP ',j,' for block ',blckn,' at ',Sys.time()),sep='\n')}
    set(test,which(test[[j]]==2),snplist[j],code.homoref[j])
    set(test,which(test[[j]]==1),snplist[j],code.het[j])
    set(test,which(test[[j]]==0),snplist[j],code.homoalt[j])
   }
#---calculate sum based on category (conserved or not conserved)
genosum.all <- test[,list(rowSums(.SD)),.SDcols=snplist]
blocklist[[blckn]] <- genosum.all
}
cat(paste0('Loop for block-wise analysis finished at ',Sys.time()),sep='\n')
blckdt <- do.call(cbind,blocklist)
fedt <- data.frame(IID=sub("\\_.*", "",rownames(gmt)),genosum=rowSums(blckdt))
write.table(fedt,paste0(outputpath,'/',outpref,'.fe.txt'),row.names=F,quote=F,sep='\t')
write.table(ncol(gmt),paste0(outputpath,'/',outpref,'.fe.nSNPs'),row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('All analysis finished at ',Sys.time()),sep='\n')
cat(paste0('Results saved to ',outputpath,'/',outpref,'.fe.txt and ',outputpath,'/',outpref,'.fe.nSNPs'),sep='\n')
