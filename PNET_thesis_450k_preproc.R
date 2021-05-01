##

################################################
#### Read in controls (Syed_GSE143209)
################################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# read in phenotype data
pd=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Syed_GSE143209_control_DNAm_data/clinical_data.csv')

# read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Syed_GSE143209_control_DNAm_data/IDAT_files")

# organize the pheno data
colnames(RGset)=paste(ss(colnames(RGset),'_',2),ss(colnames(RGset),'_',3),sep='_')
length(colnames(RGset)) == length(unique(colnames(RGset))) #TRUE
mm=match(colnames(RGset),pd$Identifier)
which(is.na(mm)) # none
pd=pd[mm,]

pd$Sex=replace(pd$Sex, which(pd$Sex=='F'),'Female')
pd$Sex=replace(pd$Sex, which(pd$Sex=='M'),'Male')

pd=data.frame(my_unique_ID=paste0('Control_', 1:nrow(pd)), age=pd$Age, sex=pd$Sex, Dx='Control', tumor_in_specimen=NA, Loc_or_Met=NA,
 Sentrix_ID=ss(pd$Identifier,'_',1), Sentrix_Position=ss(pd$Identifier,'_',2))

pData(RGset)=as(pd, "DataFrame")
colnames(RGset)=pData(RGset)$my_unique_ID

write.csv(pData(RGset),file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Syed_GSE143209_control_DNAm_data/pheno_data_for_paper.csv",row.names=FALSE)

save(RGset, file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_Syed_controls.rda")
#load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_Syed_controls.rda") # annotation: ilmn12.hg19; 64 samples

################################################
## Read in PNET data (Chan_GSE117852)
################################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_GSE117852_PNET_DNAm_data/IDAT_files")

# read in phenotype data
metadata=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_GSE117852_PNET_DNAm_data/Metadata.csv')
GSM_to_paper_IDs=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_GSE117852_PNET_DNAm_data/GSM_to_paper_IDs.csv')
	GSM_to_paper_IDs$paper_ID=toupper(ss(GSM_to_paper_IDs$paper_ID,'_methylation',1))
sentrix_info=data.frame(GSM=ss(colnames(RGset),"_",1), Sentrix_ID=ss(colnames(RGset),"_",2), Sentrix_Position=ss(colnames(RGset),"_",3))

length(ss(colnames(RGset),"_",1))==length(unique(ss(colnames(RGset),"_",1))) # TRUE
mm=match(ss(colnames(RGset),"_",1),GSM_to_paper_IDs$GSM_ID)
which(is.na(mm)) # none
GSM_to_paper_IDs=GSM_to_paper_IDs[mm,]

mm=match(GSM_to_paper_IDs$paper_ID,metadata$Paper_ID)
metadata=metadata[mm,]

metadata$Gender=replace(metadata$Gender, which(metadata$Gender=='F'),'Female')
metadata$Gender=replace(metadata$Gender, which(metadata$Gender=='M'),'Male')

pd=data.frame(my_unique_ID=metadata$Paper_ID, age=metadata$Age, sex=metadata$Gender, Dx='PNET', tumor_in_specimen=metadata$tumor_in_specimen, Loc_or_Met=metadata$Loc_or_Met,
 Sentrix_ID=sentrix_info$Sentrix_ID, Sentrix_Position=sentrix_info$Sentrix_Position)

pData(RGset)=as(pd, "DataFrame")
colnames(RGset)=pData(RGset)$my_unique_ID

# keep only tumors harvested from primary site (pancreas) (so drop tumor samples collected from metastatic sites or lymph nodes)
drop=which(pData(RGset)$tumor_in_specimen != 'Pancreas')
RGset=RGset[,-drop]

write.csv(pData(RGset),file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_GSE117852_PNET_DNAm_data/pheno_data_for_paper.csv", row.names=FALSE)

save(RGset, file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_Chan_controls.rda")
#load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_Chan_controls.rda") # annotation: ilmn12.hg19; 14 samples

################################################
## Downstream analysis
################################################

rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(calibrate)

load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_Syed_controls.rda") # annotation: ilmn12.hg19; 64 samples
RGset_Syed_controls=RGset
load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_Chan_controls.rda") # annotation: ilmn12.hg19; 14 samples
RGset_Chan_PNET=RGset

# combine the datsets

RGset=combineArrays(RGset_Syed_controls, RGset_Chan_PNET, outType = c("IlluminaHumanMethylation450k"))

write.csv(pData(RGset),file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/pheno_data_for_paper.csv", row.names=FALSE)

length(colnames(RGset)) == length(unique(colnames(RGset))) #TRUE

save(RGset, file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_combined.rda")
#load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_combined.rda")

#generate a MethylSet
MSet = preprocessRaw(RGset)

# QC
qc = getQC(MSet)
qc.df=as.data.frame(qc)

pd=pData(MSet)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/QC_mMed_vs_uMed.pdf')

plotQC(qc)

plot(qc.df$mMed, qc.df$uMed, pch=21, bg='gray', col='black', xlab='mMed', ylab='uMed', cex=1.4)

plot(qc.df$mMed, qc.df$uMed, pch=21, bg='gray', col='black',xlab='mMed', ylab='uMed', cex=1.3)
textxy(qc.df$mMed,qc.df$uMed,pd$my_unique_ID,cex =.7, offset = .7)

dev.off()

# look at Beta density plot for each sample

beta=getBeta(MSet)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/QC_density_plots.pdf')

densityPlot(MSet)

for (i in 1:ncol(beta)) {

	b=as.vector(beta[,i])
	
	if (length(which(is.na(b)) != 0)){ b=b[-which(is.na(b))] }

	plot(density(b),main=colnames(beta)[i],xlab='Beta')


}

dev.off()

# all samples look good

# check that the reported sex is correct
GMSet=mapToGenome(MSet)
check_sex=getSex(GMSet, cutoff=-2)
ps=as.vector(check_sex$predictedSex)
	ps=replace(ps, which(ps=='F'),'Female')
	ps=replace(ps, which(ps=='M'),'Male')
length(which((as.vector(pd$sex)==ps)=='FALSE')) # 0; in all samples reported sex agrees with predicted sex

GMSet=addSex(GMSet)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/sex_plots.pdf')

par(mar=c(5.1, 4.5, 4.1, 2.1))

sex=factor(pd$sex, levels=c('Female','Male'))
col=c('pink','blue')
palette(col)
plot(log2(check_sex$xMed), log2(check_sex$yMed), pch=21, bg=sex, cex=1.5, col="black", xlab=bquote(log[2](xMed)), ylab=bquote(log[2](yMed)))
legend(x='bottomleft',legend=levels(sex), col=col, pch=15, cex=1, title="Reported sex")

dev.off()

## normalize via FUNCTIONAL normalization

load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/RGset_combined.rda")
GRset=preprocessFunnorm(RGset)

pData(GRset)=as(pd, "DataFrame")
colnames(GRset)=pData(GRset)$Our_ID

save(GRset, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm.rda') # GRset

# remove probes that are common SNPs
GRset # 485512 probes; 78 samples

GRset=dropLociWithSnps_mine(object=GRset, maf=0, snpAnno = "SNPs.147CommonSingle") # 468392 probes (17120 probes dropped)

# drop sex chr
anno=getAnnotation(GRset)
GRset=GRset[! anno$chr %in% c("chrX","chrY"),] # 456928 probes

save(GRset, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

































