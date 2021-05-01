#
###########################################################################################
##### Validation analysis [case vs control] (Chan PNET dataset vs. Fadista control pancreatic islets)
###########################################################################################

###### Read in the data and peform differentrial expression analysis

# Salmon quantification files:
# /athena/masonlab/scratch/users/nai2008/PNET/Chan_etal_NatComm_2018/Salmon_quants_withGCbiasCorrection/
# /athena/masonlab/scratch/users/nai2008/Fadista_etal_control_pancreatic_islets/Salmon_quants_withGCbiasCorrection/

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
library(tximeta)

## Chan et al PNET data
salmonFiles_Chan = list.files("/athena/masonlab/scratch/users/nai2008/PNET/Chan_etal_NatComm_2018/Salmon_quants_withGCbiasCorrection", "quant.sf", recursive=T, full.names=T)
names(salmonFiles_Chan) = ss(salmonFiles_Chan,'/',10)
metadata_Chan=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_et_al_GSE118014_PNET_RNAseq_data/Metadata.csv')
SRA_and_GSM_IDs=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_et_al_GSE118014_PNET_RNAseq_data/SRA_and_GSM_IDs.csv')
GSM_to_paper_IDs=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_et_al_GSE118014_PNET_RNAseq_data/GSM_to_paper_IDs.csv')
GSM_to_paper_IDs$paper_ID=toupper(ss(GSM_to_paper_IDs$paper_ID, "\\s+", 1))

length(metadata_Chan$Paper_ID) == length(unique(metadata_Chan$Paper_ID)) #TRUE
mm=match(metadata_Chan$Paper_ID, GSM_to_paper_IDs$paper_ID)
which(is.na(mm)) # none
length(mm) == length(unique(mm)) #TRUE
metadata_Chan$GSM_ID=GSM_to_paper_IDs$GSM_ID[mm]
mm=match(metadata_Chan$GSM_ID,SRA_and_GSM_IDs$Sample.Name)
which(is.na(mm)) # none
metadata_Chan$SRA_ID=SRA_and_GSM_IDs$Run[mm]
mm=match(metadata_Chan$SRA_ID,names(salmonFiles_Chan))
which(is.na(mm)) # none
metadata_Chan$files=as.vector(salmonFiles_Chan)[mm]

metadata_Chan$Gender=replace(metadata_Chan$Gender,which(metadata_Chan$Gender=='F'),'Female')
metadata_Chan$Gender=replace(metadata_Chan$Gender,which(metadata_Chan$Gender=='M'),'Male')

# take only PNETs from pacreas (i.e. remove tumor harvested from metastatic sites)
metadata_Chan=metadata_Chan[which(metadata_Chan$tumor_in_specimen=='Pancreas'),]

metadata_Chan=data.frame(files=metadata_Chan$files , names=metadata_Chan$SRA_ID, Dx='PNET', Loc_or_Met=metadata_Chan$Loc_or_Met, Sex=metadata_Chan$Gender, Age=metadata_Chan$Age)

# Fadista control pancreatic islet data
salmonFiles_Fadista = list.files('/athena/masonlab/scratch/users/nai2008/PNET/Fadista_etal_control_pancreatic_islets/Salmon_quants_withGCbiasCorrection', "quant.sf", recursive=T, full.names=T)
names(salmonFiles_Fadista) = ss(salmonFiles_Fadista,'/',10)
metadata_Fadista_raw=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Fadista_GSE50244_control_RNAseq_data/SraRunTable.csv')
metadata_Fadista=data.frame(Sample_ID=metadata_Fadista_raw$Run, Dx='Control', Loc_or_Met=NA, Sex=simpleCap(metadata_Fadista_raw$gender), Age=metadata_Fadista_raw$Age)
mm=match(names(salmonFiles_Fadista), metadata_Fadista$Sample_ID)
metadata_Fadista=metadata_Fadista[mm,]
metadata_Fadista=cbind(as.vector(salmonFiles_Fadista),metadata_Fadista)
colnames(metadata_Fadista) = c('files','names','Dx','Loc_or_Met','Sex','Age')

metadata=rbind(metadata_Chan,metadata_Fadista)
all(ss(metadata$file,'/',10)==metadata$names) #TRUE

length(metadata$names)==length(unique(metadata$names)) #TRUE

rownames(metadata) = metadata$names
metadata$Dx=factor(metadata$Dx,levels=c('Control','PNET'))

# import data
se = tximeta(coldata=metadata, type = "salmon")
# found matching transcriptome:
# [ Ensembl - Homo sapiens - release 97 ]

# summarize transcript-level quantifications to gene-level
gse = summarizeToGene(se)

# get TPM matrix
tpm = assays(gse)$abundance

# make DESeqDataSet object
dds = DESeqDataSet(gse, design = ~ Dx)

# run SVA
library(sva)

dds <- estimateSizeFactors(dds) # using 'avgTxLength' from assays(dds), correcting for library size
dat <- counts(dds, normalized=TRUE)
idx = rowMeans(dat) > 1
dat = dat[idx, ]
mod = model.matrix(~ Dx, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svaobj = svaseq(dat, mod, mod0)
# Number of significant surrogate variables is:  20

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_SV_plots.pdf')

for (i in 2:ncol(svaobj$sv)){
	plot(svaobj$sv[,1],svaobj$sv[,i], col='black', bg=ifelse(metadata$Dx == 'PNET', "red", ifelse(metadata$Dx == 'Control',"blue","black")), pch=21, cex=1.5, xlab='SV1',ylab=paste0('SV',i))
	legend('topleft',legend=c('PNET','Control'), pch=15, col=c('red','blue'))
}

dev.off()

colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

colData(dds) = as(cbind(as.data.frame(colData(gse)),svaobj$sv),'DataFrame')

design(dds) = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + SV_10 + SV_11 + SV_12 + SV_13 + SV_14 + SV_15 + SV_16 + SV_17 + SV_18 + SV_19 + SV_20 + Dx

# run DE analysis
#dds=DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=1000)

# 236 rows did not converge in beta
# omit rows that did not converge in beta (these are typically genes with very small counts and little power)
# see https://support.bioconductor.org/p/65091/
ddsClean <- dds[which(mcols(dds)$betaConv),]

# extract results
rr=results(ddsClean, alpha=0.1, contrast=c('Dx','PNET','Control'))
# contrast = c( the name of a factor in the design formula, name of the numerator level for the fold change, name of the denominator level for the fold change)
summary(rr)
# out of 32958 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 10588, 32%
# LFC < 0 (down)     : 7290, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 3834, 12%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_Cooks_distances.pdf')

par(mar=c(8,3,2,.5))
boxplot(log10(assays(ddsClean)[["cooks"]]), range=0, xaxt='n')
axis(1, las=2, cex.axis=.3, at=1:ncol(log10(assays(ddsClean)[["cooks"]])), labels=colnames(log10(assays(ddsClean)[["cooks"]])))

dev.off()

# add gene symbols, chr, & Ensembl gene IDs
library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

mm=match(rownames(rr), genes$gene_id)
length(which(is.na(mm))) # 0
rr$chr=as.vector(genes$seqnames[mm])
rr$Ensembl=as.vector(rownames(rr))
rr$gene=as.vector(genes$gene_name[mm])

# make transformed count data, using variance stabilizing transformation (VST)
vsd = vst(ddsClean, blind=FALSE)

# save releveant data
save(se, gse, tpm, ddsClean, rr, vsd, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd

################ Downstream analysis

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
rr=rr[-which(is.na(rr$padj)),]

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

### heatmap of known M1 & M2 genes
library(pheatmap)
library(RColorBrewer)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

vst_counts = assay(vsd)

mm=match(rownames(vst_counts),genes$gene_id)
which(is.na(mm)) # none
rownames(vst_counts) = genes$gene_name[mm]

mm_M1=match(M1_marker_panel$gene_symbol, rownames(vst_counts))
which(is.na(mm)) # none
mm_M2=match(M2_marker_panel$gene_symbol, rownames(vst_counts))
which(is.na(mm)) # none

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_M1_and_M2_heatmaps.pdf')

#annotation_col=data.frame(Dx=colData(ddsClean)[,'Dx'])
annotation_col=data.frame(Dx=colData(ddsClean)[,'Dx'],Met.Status=colData(ddsClean)[,'Loc_or_Met'])
rownames(annotation_col)=as.vector(colData(ddsClean)[,'names'])

pheatmap(vst_counts[mm_M1,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='M1 marker panel (VST transformed counts)', 
	cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=3)

pheatmap(vst_counts[mm_M1,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='M1 marker panel (VST transformed counts)', 
	cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, fontsize_col=3)

pheatmap(vst_counts[mm_M2,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='M2 marker panel (VST transformed counts)', 
	cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, scale='row', fontsize_col=3, cutree_cols=5)

pheatmap(vst_counts[mm_M2,], color=colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(3000), main='M2 marker panel (VST transformed counts)', 
	cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=annotation_col, fontsize_col=3)

dev.off()

### make table of DE genes

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=sig_results$gene, chr=sig_results$chr, Ensembl=sig_results$Ensembl, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

nrow(out) # 17878

# How many DE genes are TFs?

out$TF=FALSE

TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)

mm=match(out$Ensembl,TFs$Ensembl_ID)
out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) # 2153

save(out, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DE_genes_byDx.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DE_genes_byDx.rda')

# print table of DE genes
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Chan_and_Fadista_datasets_DE_genes_byDx.csv", row.names=FALSE)

### calculate concordance scores with the mouse transcriptome for M1 and M2 genes

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_M1_and_M2_concordance_scores.pdf')

par(mar=c(5.1, 4.7, 4.1, 2.1))

# M1
mm = match(genes_upregulated_in_M1_and_downregulated_in_M2$Ensembl, rr$Ensembl)
which(is.na(mm)) # none
rr_M1 = rr[mm,]
disco.score = rr_M1$log2FoldChange * log2(genes_upregulated_in_M1_and_downregulated_in_M2$FC_M1vsM0) * abs(log10(rr_M1$padj) + log10(genes_upregulated_in_M1_and_downregulated_in_M2$p.value_M1vsM0))
oo=order(disco.score)

length(which(disco.score>0))/length(disco.score) # 0.6097561

plot(disco.score[oo],1:length(disco.score[oo]),type='n', axes=FALSE, ylab=NA, xlab='Concordance score', main='Concordance between M1 genes')
axis(1)
par(las=2)
axis(2, at=1:length(disco.score[oo]),labels=rr_M1$gene,cex.axis=0.7)
box(lty='solid')
for(i in 1:length(disco.score[oo])){abline(i,0, col = "lightgray", lty = "dotted")}
points(disco.score[oo],1:length(disco.score[oo]), pch=23, col='black', bg='grey')
abline(v=0,col='black',lty='solid',lwd=0.7)

# M2
mm = match(genes_upregulated_in_M2_and_downregulated_in_M1$Ensembl, rr$Ensembl)
which(is.na(mm)) # none
rr_M2 = rr[mm,]
disco.score = rr_M2$log2FoldChange * log2(genes_upregulated_in_M2_and_downregulated_in_M1$FC_M2vsM0) * abs(log10(rr_M2$padj) + log10(genes_upregulated_in_M2_and_downregulated_in_M1$p.value_M2vsM0))
oo=order(disco.score)

length(which(disco.score>0))/length(disco.score) # 0.4545455

par(las=1)
plot(disco.score[oo],1:length(disco.score[oo]),type='n', axes=FALSE, ylab=NA, xlab='Concordance score', main='Concordance between M2 genes')
axis(1)
par(las=2)
axis(2, at=1:length(disco.score[oo]),labels=rr_M2$gene,cex.axis=0.7)
box(lty='solid')
for(i in 1:length(disco.score[oo])){abline(i,0, col = "lightgray", lty = "dotted")}
points(disco.score[oo],1:length(disco.score[oo]), pch=23, col='black', bg='grey')
abline(v=0,col='black',lty='solid',lwd=0.7)

dev.off()

### Perform gene set over-representation analysis (ORA)

library(goseq)
load('/athena/masonlab/scratch/users/nai2008/items_for_goseq_analysis.rda') # gene2cat_GOandKEGG, KEGG_term_names, median_tx_lengths, cat2gene_GO, cat2gene_KEGG

indicator=rep(0, times=nrow(rr))
indicator[which(rr$padj<=0.1)]=1
aa=indicator
names(aa)=rr$Ensembl

mm = match(names(aa), median_tx_lengths$gene_EnsemblID)
bias.data = median_tx_lengths$median_length[mm]
pwf = nullp(aa, 'hg38', 'ensGene', bias.data = bias.data)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_goseq_pwf_plot.pdf')
plotPWF(pwf)
dev.off()

GO.KEGG.wall=goseq(pwf,"hg38","ensGene", gene2cat = gene2cat_GOandKEGG, test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"))
GO.KEGG.wall$over_represented_FDR=p.adjust(GO.KEGG.wall$over_represented_pvalue, method="BH")

GO.KEGG.wall$ontology[grep('path:hsa', GO.KEGG.wall$category)]='KEGG'

index = grep('path:hsa', GO.KEGG.wall$category)

for (i in 1:length(index)){
	mm=match(GO.KEGG.wall$category[index[i]], KEGG_term_names$KEGG_ID)
	GO.KEGG.wall$term[index[i]] = KEGG_term_names$KEGG_term[mm]
}

length(which(GO.KEGG.wall$over_represented_FDR<=0.1)) #115
GO.KEGG.wall_sig = GO.KEGG.wall[which(GO.KEGG.wall$over_represented_FDR<=0.1),]

# # Add DE genes in each GO/KEGG category

# GO.KEGG.wall_sig_withoutGenes = GO.KEGG.wall_sig

# library(EnsDb.Hsapiens.v86)
# edb = EnsDb.Hsapiens.v86
# ens.gene.map = genes(edb, columns = c("gene_id", "gene_name"), return.type="data.frame")

# length(names(cat2gene_GO)) == length(unique(names(cat2gene_GO))) #TRUE
# length(names(cat2gene_KEGG)) == length(unique(names(cat2gene_KEGG))) #TRUE

# GO.KEGG.wall_sig$genes_Ensembl=NA
# GO.KEGG.wall_sig$genes=NA

# for (i in 1:nrow(GO.KEGG.wall_sig)){

# 	cat=GO.KEGG.wall_sig$category[i]

# 	if (length(grep('GO',cat)) == 1){

# 		m.cat=match(cat, names(cat2gene_GO))

# 		if(is.na(m.cat)){print('error: m.cat does not match (GO)')} else {

# 			possible_genes=cat2gene_GO[[m.cat]]

# 			m.genes=match(possible_genes,names(aa))

# 			if( length(which(!is.na(m.genes)))==0 ){print('error: m.genes are all <NA> (GO)')} else {

# 				if (length(which(is.na(m.genes)))>0){ possible_genes= possible_genes[-which(is.na(m.genes))] }

# 				m.genes=match(possible_genes,names(aa))
# 				subset=aa[m.genes]
# 				DE_genes=subset[which(subset==1)]
# 				GO.KEGG.wall_sig$genes_Ensembl[i]=paste(names(DE_genes),collapse='')

# 				m.ens=match(names(DE_genes),ens.gene.map$gene_id)
# 				GO.KEGG.wall_sig$genes[i]=paste(ens.gene.map$gene_name[m.ens],collapse='')
# 				}
# 		}
# 	} else if (length(grep('path:hsa',cat)) == 1){

# 		m.cat=match(cat, names(cat2gene_KEGG))

# 		if(is.na(m.cat)){print('error: m.cat does not match (KEGG)')} else {

# 			possible_genes=cat2gene_KEGG[[m.cat]]

# 			m.genes=match(possible_genes,names(aa))

# 			if(length(which(!is.na(m.genes) == 0))){print('error: m.genes are all <NA> (KEGG)')} else{

# 				if (length(which(is.na(m.genes)))>0){ possible_genes= possible_genes[-which(is.na(m.genes))] }

# 				m.genes=match(possible_genes,names(aa))
# 				subset=aa[m.genes]
# 				DE_genes=subset[which(subset==1)]
# 				GO.KEGG.wall_sig$genes_Ensembl[i]=paste(names(DE_genes),collapse=';')

# 				m.ens=match(names(DE_genes),ens.gene.map$gene_id)
# 				GO.KEGG.wall_sig$genes[i]=paste(ens.gene.map$gene_name[m.ens],collapse=';')
# 			}
# 		}
# 	}
# }

write.csv(GO.KEGG.wall_sig, file = '/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Chan_and_Fadista_datasets_DEbyDx_ORA.csv')

# make ORA plots

ora_to_plot = GO.KEGG.wall_sig[which(GO.KEGG.wall_sig$ontology %in% c('BP','KEGG')),]

if(length(which(ora_to_plot$ontology == 'KEGG')) > 0 ){
	ora_to_plot$term[which(ora_to_plot$ontology == 'KEGG')] = paste0(ora_to_plot$term[which(ora_to_plot$ontology == 'KEGG')], ' (KEGG)')
}

ora_to_plot = ora_to_plot[order(ora_to_plot$over_represented_FDR, decreasing = TRUE),]

ora_to_plot$neg_log_over_represented_FDR = -log(ora_to_plot$over_represented_FDR)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_ORA_plot.pdf')
par(mar=c(5.1, 15, 4.1, 0))
barplot(ora_to_plot$neg_log_over_represented_FDR, names.arg= ora_to_plot$term, xlim=c(0,5+round(max(ora_to_plot$neg_log_over_represented_FDR))), horiz=TRUE, col='black', xlab = '-Log(FDR)', ylab=NA, las=1, cex.names=0.5)
dev.off()

### make TPM plots of relevant genes

genes_to_plot=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/genes_of_interest.csv', header=TRUE)

mm=match(genes_to_plot$ensembl, rr$Ensembl)
if(length(which(is.na(mm))) != 0) { genes_to_plot=genes_to_plot[-which(is.na(mm)),] }
mm=match(genes_to_plot$ensembl, rr$Ensembl)
genes_of_interest_rr_subset = rr[mm,]
genes_of_interest_rr_subset$info = genes_to_plot$info
save(genes_of_interest_rr_subset, file = '/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_genes_of_interest_rr_subset.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_genes_of_interest_rr_subset.rda') # genes_of_interest_rr_subset

nrow(genes_to_plot) == length(unique(genes_to_plot$ensembl)) #TRUE

mm = match(genes_to_plot$ensembl, genes$gene_id)
which(is.na(mm)) # none
genes_to_plot$gene_symbol=genes$gene_name[mm]

log2_tpm_plus_one=log2(tpm+1)

all(colnames(log2_tpm_plus_one)==as.vector(ddsClean$names)) #TRUE

mycol=as.vector(ddsClean$Dx)
mycol[which(mycol=="Control")]='blue'
mycol[which(mycol=="PNET")]='red'

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_genes_of_interest_TPMs.pdf')

for (i in 1: length(genes_to_plot$ensembl)){

	par(mar=c(5.1,5.3,4.1,2.1))

	index=which(rr$Ensembl==genes_to_plot$ensembl[i])

	zag1=paste0(rr$gene[index],' (', genes_to_plot$info[i], ')')
	zag2=as.expression(bquote(log[2]~"FC" == .(signif(rr$log2FoldChange[index],2))))
	zag3=paste0("FDR = ",signif(rr$padj[index],2))	

	log2_tpm_plus1_subset=as.vector(log2_tpm_plus_one[which(rownames(log2_tpm_plus_one)==genes_to_plot$ensembl[i]),])
	x=factor(as.vector(ddsClean$Dx),levels=c('Control','PNET'))

	boxplot(as.vector(log2_tpm_plus1_subset)~x, , xlab='Disease state', ylab= as.expression(bquote(log[2]~"(TPM+1)")), main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(as.vector(log2_tpm_plus1_subset) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topleft',legend=c(zag2,zag3), bty='n')

}

dev.off()

### volcano plot

genes_to_plot=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/genes_of_interest.csv', header=TRUE)

mm = match(genes_to_plot$ensembl, genes$gene_id)
which(is.na(mm)) # none
genes_to_plot$gene_symbol=genes$gene_name[mm]

genes.to.show=as.vector(genes_to_plot[which(genes_to_plot$info %in% c('M1 marker','M2 marker')),])
genes.to.show$col=NA
genes.to.show$col[which(genes_to_plot$info == 'M1 marker')]='blue'
genes.to.show$col[which(genes_to_plot$info == 'M2 marker')]='green'

mm=match(genes.to.show$gene_symbol,rr$gene)
genes.to.show=genes.to.show[-which(is.na(mm)),]

library(EnhancedVolcano)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_volcano_plot_DEbyDx.pdf')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"FC"~"|"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"FC"~"|"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano(toptable = as.data.frame(rr[,c(2,6)]),
	lab = rr$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='Validation analysis',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	#xlim=c(-35,41),
	#ylim=c(-5.5,305),
	shape=19,
	legendLabSize=8,
	legendIconSize=15
	)

EnhancedVolcano(toptable=as.data.frame(rr[,c(2,6)]),
	lab = rr$gene,
	x = 'log2FoldChange',
	y = 'padj',
	selectLab=genes.to.show$gene_symbol,
	labCol='black',
	labFace='bold',
	#boxedLabels=TRUE,
	labSize = 3,
	drawConnectors=TRUE,
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M1 and M2 marker genes labeled',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	xlim=c(-30,30),
	ylim=c(-2,75),
	shape=19,
	legendLabSize=8,
	legendIconSize=15
	)

EnhancedVolcano(toptable=as.data.frame(rr[,c(2,6)]),
	lab = rr$gene,
	x = 'log2FoldChange',
	y = 'padj',
	selectLab=genes.to.show$gene_symbol[grep('M1',genes.to.show$info)],
	labCol='black',
	labFace='bold',
	#boxedLabels=TRUE,
	labSize = 3,
	drawConnectors=TRUE,
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M1 marker genes lableled',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	xlim=c(-30,30),
	ylim=c(-2,75),
	shape=19,
	#labhjust = 0,
	#labvjust = -1,
	legendLabSize=8,
	legendIconSize=15
	)

EnhancedVolcano(toptable=as.data.frame(rr[,c(2,6)]),
	lab = rr$gene,
	x = 'log2FoldChange',
	y = 'padj',
	selectLab=genes.to.show$gene_symbol[grep('M2',genes.to.show$info)],
	labCol='black',
	labFace='bold',
	#boxedLabels=TRUE,
	labSize = 3,
	drawConnectors=TRUE,
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M2 marker genes lableled',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="none",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1/2,
	xlim=c(-20,20),
	ylim=c(-2,50),
	shape=19,
	#labhjust=-3,
	#labvjust = -2,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"Fold"~"Change"~"|"~">"~"1"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"Fold"~"Change"~"|"~">"~"1"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=19, col=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()

### in silico immune cell sorting

library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

# make TPM matrix
mm=match(rownames(tpm),genes$gene_id)
which(is.na(mm)) # none
rownames(tpm)=genes$gene_name[mm]

pnet_tpm_mat=tpm[,which(as.vector(ddsClean$Dx)=='PNET')]
control_tpm_mat=tpm[,which(as.vector(ddsClean$Dx)=='Control')]

# make TPM matrix with removed duplicated genes
tpm_dup_rm=tpm[-c(which(duplicated(rownames(tpm))==TRUE),which(duplicated(rownames(tpm),fromLast=TRUE)==TRUE)),]
pnet_tpm_dup_rm = tpm_dup_rm[,which(as.vector(ddsClean$Dx)=='PNET')]
control_tpm_dup_rm = tpm_dup_rm[,which(as.vector(ddsClean$Dx)=='Control')]

## quanTIseq
res_quantiseq_pnet = deconvolute(pnet_tpm_mat, method="quantiseq", tumor = TRUE)
res_quantiseq_control = deconvolute(control_tpm_mat, method="quantiseq", tumor = FALSE)

## EPIC
res_epic_pnet = deconvolute(pnet_tpm_dup_rm, method="epic", tumor = TRUE)
# Warning messages:
# 1: In (function (bulk, reference = NULL, mRNA_cell = NULL, mRNA_cell_sub = NULL,  :
#   The optimization didn't fully converge for some samples:
# SRR7634542
pnet_tpm_dup_rm_excluding_SRR7634542=pnet_tpm_dup_rm[,-which(colnames(pnet_tpm_dup_rm)=='SRR7634542')]

res_epic_pnet = deconvolute(pnet_tpm_dup_rm_excluding_SRR7634542, method="epic", tumor = TRUE)
res_epic_control = deconvolute(control_tpm_dup_rm, method="epic", tumor = FALSE)

## CIBERSORT (abs. mode)
set_cibersort_binary("/athena/masonlab/scratch/users/nai2008/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("/athena/masonlab/scratch/users/nai2008/CIBERSORT/LM22.txt")

res_cibersort = deconvolute(tpm_dup_rm, method="cibersort_abs")

## MCP-counter
res_mcpCounter = deconvolute(tpm, method="mcp_counter")

## xCell
res_xCell = deconvolute(tpm, method='xcell')

## make BAR plots comparing cell type fractions

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_inSilico_deconvolution.pdf')

ggplot = function(...) { 	ggplot2::ggplot(...) + theme_bw() }

#quantiseq

res_quantiseq_pnet %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq_pnet))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (quanTIseq); PNET Samples')

res_quantiseq_pnet %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq_pnet))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (quanTIseq); PNET Samples') +
    theme(legend.position = "none")

res_quantiseq_control %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq_control))) +
    ylab('Cell fraction') +
    xlab('Sample')+
    theme(axis.text.y = element_text(size = rel(0.6))) +
    guides(fill=guide_legend(title="Cell Types ")) +
    ggtitle('Immune cell type deconvolution (quanTIseq); Control pancreatic islets')


res_quantiseq_control %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq_control))) +
    ylab('Cell fraction') +
    xlab('Sample')+
    theme(axis.text.y = element_text(size = rel(0.6))) +
    guides(fill=guide_legend(title="Cell Types ")) +
    ggtitle('Immune cell type deconvolution (quanTIseq); Control pancreatic islets') +
    theme(legend.position = "none")

#epic

res_epic_pnet %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_epic_pnet))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); PNET Samples')

res_epic_pnet %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_epic_pnet))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); PNET Samples') +
    theme(legend.position = "none")

res_epic_control %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_epic_control))) +
    ylab('Cell fraction') +
    xlab('Sample')+
    theme(axis.text.y = element_text(size = rel(0.6))) +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); Control pancreatic islets')

res_epic_control %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_epic_control))) +
    ylab('Cell fraction') +
    xlab('Sample')+
    theme(axis.text.y = element_text(size = rel(0.6))) +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); Control pancreatic islets') +
    theme(legend.position = "none")

dev.off()

## make BOX plots comparing cell type abundance

# quantiseq
res_quantiseq_pnet_df=as.data.frame(res_quantiseq_pnet)
res_quantiseq_control_df=as.data.frame(res_quantiseq_control)

all(res_quantiseq_pnet_df$cell_type==res_quantiseq_control_df$cell_type) #TRUE

p.wilcox.test=vector()
for (i in 1:nrow(res_quantiseq_pnet_df)){
	p=wilcox.test(as.numeric(res_quantiseq_control_df[i,2:ncol(res_quantiseq_control_df)]), as.numeric(res_quantiseq_pnet_df[i,2:ncol(res_quantiseq_pnet_df)]))
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')
	
pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_deconvolution_cellType_comparisons_quantiseq.pdf')

for (i in 1:nrow(res_quantiseq_pnet_df)){

	cell_type_data=c(as.numeric(res_quantiseq_pnet_df[i,2:ncol(res_quantiseq_pnet_df)]),as.numeric(res_quantiseq_control_df[i,2:ncol(res_quantiseq_control_df)]))
	dx=c(rep('PNET',times=ncol(res_quantiseq_pnet_df)-1),rep('Control',times=ncol(res_quantiseq_control_df)-1))
	dx=factor(dx,levels=c('Control','PNET'))

	median_PNET=median(as.numeric(res_quantiseq_pnet_df[i,2:ncol(res_quantiseq_pnet_df)]))
	median_control=median(as.numeric(res_quantiseq_control_df[i,2:ncol(res_quantiseq_control_df)]))

	mycol=as.vector(dx)
	mycol[which(mycol=="Control")]='blue'
	mycol[which(mycol=="PNET")]='red'
 
	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=res_quantiseq_pnet_df$cell_type[i]
	zag2=paste0('{ Control median = ', signif(median_control,2), '; PNET median = ', signif(median_PNET,2),' }')

	boxplot(cell_type_data~dx, , xlab='Disease state', ylab= 'Cell fraction', main=c(zag1,zag2), cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(dx), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topleft',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#epic
res_epic_pnet_df = as.data.frame(res_epic_pnet)
res_epic_control_df = as.data.frame(res_epic_control)

cct = intersect(res_epic_control_df$cell_type, res_epic_pnet_df$cell_type)

res_epic_pnet_df = res_epic_pnet_df[match(cct,res_epic_pnet_df$cell_type),]
res_epic_control_df = res_epic_control_df[match(cct,res_epic_control_df$cell_type),]

all(res_epic_pnet_df$cell_type==res_epic_control_df$cell_type) #TRUE

p.wilcox.test=vector()
for (i in 1:nrow(res_epic_pnet_df)){
	p=wilcox.test(as.numeric(res_epic_control_df[i,2:ncol(res_epic_control_df)]), as.numeric(res_epic_pnet_df[i,2:ncol(res_epic_pnet_df)]))
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_deconvolution_cellType_comparisons_epic.pdf')

for (i in 1:nrow(res_epic_pnet_df)){

	cell_type_data=c(as.numeric(res_epic_pnet_df[i,2:ncol(res_epic_pnet_df)]),as.numeric(res_epic_control_df[i,2:ncol(res_epic_control_df)]))
	dx=c(rep('PNET',times=ncol(res_epic_pnet_df)-1),rep('Control',times=ncol(res_epic_control_df)-1))
	dx=factor(dx,levels=c('Control','PNET'))

	median_PNET=median(as.numeric(res_epic_pnet_df[i,2:ncol(res_epic_pnet_df)]))
	median_control=median(as.numeric(res_epic_control_df[i,2:ncol(res_epic_control_df)]))

	mycol=as.vector(dx)
	mycol[which(mycol=="Control")]='blue'
	mycol[which(mycol=="PNET")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=res_epic_pnet_df$cell_type[i]
	zag2=paste0('{ Control median = ', signif(median_control,2), '; PNET median = ', signif(median_PNET,2),' }')

	boxplot(cell_type_data~dx, , xlab='Disease state', ylab= 'Cell fraction', main=c(zag1,zag2), cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(dx), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#CIBERSORT (abs. mode)
res_cibersort = as.data.frame(res_cibersort)

cell_type = as.vector(res_cibersort$cell_type)
abundance_matrix = res_cibersort[,2:ncol(res_cibersort)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Dx)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'PNET')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Control')]))
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_deconvolution_cellType_comparisons_cibersort.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	dx=factor(colnames(abundance_matrix),levels=c('Control','PNET'))

	median_PNET=median(cell_type_data[which(dx=='PNET')])
	median_control=median(cell_type_data[which(dx=='Control')])

	mycol=as.vector(dx)
	mycol[which(mycol=="Control")]='blue'
	mycol[which(mycol=="PNET")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]
	zag2=paste0('{ Control median = ', signif(median_control,2), '; PNET median = ', signif(median_PNET,2),' }')

	boxplot(cell_type_data~dx, , xlab='Disease state', ylab= 'Cell type abundance (arbitrary units)', main=c(zag1,zag2), cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(dx), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#MCP-counter
res_mcpCounter = as.data.frame(res_mcpCounter)

cell_type = as.vector(res_mcpCounter$cell_type)
abundance_matrix = res_mcpCounter[,2:ncol(res_mcpCounter)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Dx)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'PNET')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Control')]))
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_deconvolution_cellType_comparisons_mcpCounter.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	dx=factor(colnames(abundance_matrix),levels=c('Control','PNET'))

	median_PNET=median(cell_type_data[which(dx=='PNET')])
	median_control=median(cell_type_data[which(dx=='Control')])

	mycol=as.vector(dx)
	mycol[which(mycol=="Control")]='blue'
	mycol[which(mycol=="PNET")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]
	zag2=paste0('{ Control median = ', signif(median_control,2), '; PNET median = ', signif(median_PNET,2),' }')

	boxplot(cell_type_data~dx, , xlab='Disease state', ylab= 'Cell population abundance (arbitrary units)', main=c(zag1,zag2), cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(dx), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#xCell
res_xCell = as.data.frame(res_xCell)

cell_type = as.vector(res_xCell$cell_type)
abundance_matrix = res_xCell[,2:ncol(res_xCell)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Dx)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'PNET')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Control')]))
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Chan_and_Fadista_datasets_deconvolution_cellType_comparisons_xCell.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	dx=factor(colnames(abundance_matrix),levels=c('Control','PNET'))

	median_PNET=median(cell_type_data[which(dx=='PNET')])
	median_control=median(cell_type_data[which(dx=='Control')])

	mycol=as.vector(dx)
	mycol[which(mycol=="Control")]='blue'
	mycol[which(mycol=="PNET")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]
	zag2=paste0('{ Control median = ', signif(median_control,2), '; PNET median = ', signif(median_PNET,2),' }')

	boxplot(cell_type_data~dx, , xlab='Disease state', ylab= 'Enrichment score (arbitrary units)', main=c(zag1,zag2), cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(dx), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()



