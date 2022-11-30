
#!/usr/bin/env Rscript

# Script to score dosage bias per sample

#Inputs:
# 1) a binCov matrix subsetted to the testing bins of interest
# 2) a list of bins with signed weights per bin

####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)

## VIASH START
#DEV TEST RUN (on local machine)
# batch.ideal <- 100
WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_analysis/WGD_version2_Dec2017/"
par <- list(
  "matrix" = paste0(WRKDIR, "WGD_training.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz"),
  "scoring_mask" = paste0(WRKDIR, "WGD_scoring_mask.6F_adjusted.100bp.h37.bed"),
  "scores" = '~/scratch/WGDscore_testing/scores.txt.gz',
  "plot" = '~/scratch/WGDscore_testing/plot.pdf',
  "gzip" = TRUE,
  "outliers" = FALSE
)
meta <- list(
  "resources_dir" = 'src/wgd_score',
)
## VIASH END

# import helper functions
source(paste0(meta["resources_dir"], "/helper.R"))


#####PART 1: DATA PROCESSING#####
#Read scoring bins
WGD.bins <- read.table(par$scoring_mask,header=T,comment.char="",sep="\t")
if(ncol(WGD.bins) != 4){
  stop("WGD bins file must have four columns; see documentation.")
}
colnames(WGD.bins) <- c("chr","start","end","weight")

#Read, normalize, and clean coverage data
cov <- readMatrix(par$matrix)
cov <- normalizeContigsPerSample(cov,exclude=c("X","Y"),ploidy=2)
cov <- filterZeroBins(cov)
colnames(cov)[1:3] <- c("chr","start","end")

#Subset coverage data to only include WGD bins
cov <- merge(WGD.bins[,1:3],cov,by=1:3)

#Subset WGD bins to only include relevant coverage data
WGD.bins <- merge(WGD.bins,cov[,1:3],by=1:3)

#####PART 2: DOSAGE SCORING PER SAMPLE#####
#Compute scores and return as data frame
scores <- scoreSamples(cov,WGD.bins)

#Write scores to file
colnames(scores)[1] <- "#ID"
write.table(scores,par$scores,col.names=T,row.names=F,sep="\t",quote=F)
if(par$gzip){
  system(paste0("gzip -f ",OUTDIR,"/WGD_scores.txt"),intern=F,wait=F)
}
colnames(scores)[1] <- "ID"

#####PART 3: PLOT SCORE DISTRIBUTIONS#####
#Only run if optioned
if(!is.null(par$plot)){
  pdf(par$plot,height=4,width=8)
  plotScores(scores,label.outliers=par$outliers)
  dev.off()
}

