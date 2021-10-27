# Genome-wide-heterozygosity and fitness
R scripts to estimate heterozygosity across a set of loci or SNPs, which can be used to test associations with fitness.

Overview: the Rscripts implemenat the analysis which generates a quantitative fixed effects of genomic heterozygosity (gHet) across the genome  or a selection of SNPs, described in Xiang et al 2021 (https://www.biorxiv.org/content/10.1101/2021.04.19.440546v1). The purpose of this fixed effects is to have it fitted to a linear (mix) model to test if this heterzygosity is associated with a y, such as a complex trait. 

If the choice of SNPs used to generate gHet involves a set of highly conserved SNPs across species, we can test if a trait is related to fitness. In Xiang et al 2021 (https://www.biorxiv.org/content/10.1101/2021.04.19.440546v1), we partition the genome into SNPs at conserved sites and SNPs at the remaining sites, then we generated 2 fixed effects of gHet estimated from the 2 sets of SNPs. Then we fit these 2 fixed effects of gHet into a linear mixed model (along with other fixed and random effects) to see how the trait is associated with them. We found traits like height or fertility is only significnatly associated with gHet from conserved sites but not with gHet from the remaining sites. Therefore, traits like height or fertility are more likely to be related to fitness.

Two Rscripts are given to estimate gHet. gHet.R is to estimate gHet and we recommend doing this per chromosome if you have dense genotype data. cbn.gHet.R is to combine results generated from gHet.R per chromosome into a single fixed effects of gHet.

To run the Rscript you will need R packages of 'data.table' and 'BEDMatrix' pre-installed in your R system. Also, the script is designed to work with plink binary files. 

Example plink datasets (https://github.com/rxiangr/Genome-wide-heterozygosity/tree/main/testplink and https://github.com/rxiangr/Genome-wide-heterozygosity/tree/main/hetfeout) are provided with a tutorial of running the analysis: https://github.com/rxiangr/Genome-wide-heterozygosity/blob/tutorial/README.md 

Reference:
Xiang R, Breen E, Bolormaa S, Vander Jagt C, Chamberlain A, Macleod I, Goddard M. Mutant alleles differentially shape fitness and other complex traits in cattle. Communications Biology (In Press). 2021
Preprint: https://www.biorxiv.org/content/10.1101/2021.04.19.440546v2
