# Usage of gHet.R and cbn.gHet.R:

#==**note that to run gHet.R and cbn.gHet.R you need to install R packages of 'data.table' and 'BEDMatrix'**

#module load R/4.0.0-foss-2020a

**Rscript gHet.R <../yourpath/prefix of plink file> <../yourpath/full name of plink file.frq> <../yourpath of output> <"your outputname">**

You can run gHet.R using example data (in folder testplink). There are plink binary (*bim/fam/bed) and frequency file (*.frq) which can be calculated using plink:

**./plink --cow --keep-allele-order  --bfile testplink/test.chr25 --freq --out test.chr25**

Then you using gHet.R to read in plink genotypes and frequency data:

**Rscript gHet.R testplink/test.chr25 testplink/test.chr25.frq hetfeout test.chr25**

This will generate 2 files: **hetfeout/test.chr25.fe.txt** (also in example data file in the main folder)

IID     genosum

200001697       -4740.504

200001718       -6577.204

200004308       -4515.714

200016610       -3639.084

200020986       -4375.444


which contains two columns. the 1st column is the individual ID and 2nd column is the sum of heterozygosity across SNPs for each individual in the supplied plink file.

and **hetfeout/test.chr25.fe.nSNPs**

10000

which contains a value of number of SNPs in this analysis. Sum of heterozygosity and number of SNPs used can be used to calculate the mean of heterozygosity.

If you repeat the above process for another chromosome, say chromosome 24, you will get another set of results, e.g., **hetfeout/test.chr24.fe.txt** and **hetfeout/test.chr24.fe.nSNPs**. 

We then use the cbn.gHet.R to combine the per-chromosome results:

 **Rscript cbn.gHet.R <../yourpathOf*fe.txt>  <../yourpath of output>  <"your outputname">** 

Using the example data:

**Rscript cbn.gHet.R hetfeout . test**

This will generate an output **./test.cbfe.qc.txt**

200001697       200001697       -0.3524941

200001718       200001718       -0.3952831

200004308       200004308       -0.3207196

200016610       200016610       -0.3021041

200020986       200020986       -0.4385046

201002087       201002087       -0.3294511

The 1st two columns are repeated individual IDs and the 3rd column is the mean heterozygosity for individuals. 
This file is compatible to quantitative fixed effects format of GCTA --qcovar (https://cnsgenomics.com/software/gcta/#GREMLanalysis) to be modelled with other random and/or fixed effects. 
