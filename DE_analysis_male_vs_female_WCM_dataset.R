#
###########################################################################################
##### Discovery analysis [male vs female] (WCM PNET dataset)
###########################################################################################

###### Read in the data and peform differentrial expression analysis

# Salmon quantification files:
# /athena/masonlab/scratch/users/nai2008/PNET/WCM_PNET_RNAseq/Salmon_quants_withGCbiasCorrection/

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
library(tximeta)

# WCM PNET Data
salmonFiles_WCM = list.files("/athena/masonlab/scratch/users/nai2008/PNET/WCM_PNET_RNAseq/Salmon_quants_withGCbiasCorrection", "quant.sf", recursive=T, full.names=T)
names(salmonFiles_WCM) = ss(salmonFiles_WCM,'/',10)
metadata_WCM=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_WCM_PNET_RNAseq_data/PNET_metadata.csv')
mm=match(names(salmonFiles_WCM), metadata_WCM$Sample_ID)
metadata_WCM = metadata_WCM[mm,]
metadata_WCM = cbind(as.vector(salmonFiles_WCM),metadata_WCM)
colnames(metadata_WCM) = c('files','names','Dx','Loc_or_Met','Sex','Age')

metadata=metadata_WCM
all(ss(metadata$file,'/',10)==metadata$names) #TRUE

metadata$names[grep('Sample',metadata$names)]=ss(metadata$names[grep('Sample',metadata$names)],'Sample_',2)
metadata$names[grep('BMF',metadata$names)]=ss(metadata$names[grep('_BMF',metadata$names)],'_BMF',1)
length(metadata$names)==length(unique(metadata$names)) #TRUE

rownames(metadata) = metadata$names

metadata$Sex=replace(metadata$Sex,which(metadata$Sex=='F'),'Female')
metadata$Sex=replace(metadata$Sex,which(metadata$Sex=='M'),'Male')
metadata$Sex=factor(metadata$Sex, levels=c('Male','Female'))
metadata$Loc_or_Met=factor(metadata$Loc_or_Met,levels=c('Localized','Metastatic'))

# import data
se = tximeta(coldata=metadata, type = "salmon")
# found matching transcriptome:
# [ Ensembl - Homo sapiens - release 97 ]

# summarize transcript-level quantifications to gene-level
gse = summarizeToGene(se)

# get TPM matrix
tpm = assays(gse)$abundance

# make DESeqDataSet object
dds = DESeqDataSet(gse, design = ~ Sex)

# run SVA
library(sva)

dds <- estimateSizeFactors(dds) # using 'avgTxLength' from assays(dds), correcting for library size
dat <- counts(dds, normalized=TRUE)
idx = rowMeans(dat) > 1
dat = dat[idx, ]
mod = model.matrix(~ Sex, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svaobj = svaseq(dat, mod, mod0)
# Number of significant surrogate variables is: 6

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_SV_plots.pdf')

for (i in 2:ncol(svaobj$sv)){
	plot(svaobj$sv[,1],svaobj$sv[,i], col='black', bg=ifelse(metadata$Sex == 'Male', "blue", ifelse(metadata$Sex == 'Female',"red","black")), pch=21, cex=1.5, xlab='SV1',ylab=paste0('SV',i))
	legend('topleft',legend=c('Male','Female'), pch=15, col=c('blue','red'))
}

dev.off()

colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

colData(dds) = as(cbind(as.data.frame(colData(gse)),svaobj$sv),'DataFrame')

design(dds) = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + Sex

# run DE analysis
dds=DESeq(dds)

ddsClean = dds

# extract results
rr=results(ddsClean, alpha=0.1, contrast=c('Sex','Male','Female')) 
# contrast = c( the name of a factor in the design formula, name of the numerator level for the fold change, name of the denominator level for the fold change  )
summary(rr)
# out of 32816 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 338, 1%
# LFC < 0 (down)     : 221, 0.67%
# outliers [1]       : 0, 0%
# low counts [2]     : 1237, 3.8%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_Cooks_distances.pdf')

par(mar=c(8,5,2,2))
boxplot(log10(assays(ddsClean)[["cooks"]]), range=0, las=2, cex.axis=.8)

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
save(se, gse, tpm, ddsClean, rr, vsd, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DESeq2_DEbySex_BUNDLE.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, ddsClean, rr, vsd

################ Downstream analysis

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
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

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_M1_and_M2_heatmaps.pdf')

#annotation_col=data.frame(Dx=colData(ddsClean)[,'Dx'])
annotation_col=data.frame(Met.Status=colData(ddsClean)[,'Loc_or_Met'])
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

nrow(out) #559

# How many DE genes are TFs?

out$TF=FALSE

TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)

mm=match(out$Ensembl,TFs$Ensembl_ID)
out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) #55

save(out, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DE_genes_bySex.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DE_genes_bySex.rda')

# print table of DE genes
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Male_vs_Female_WCM_dataset_DE_genes_bySex.csv", row.names=FALSE)

### calculate concordance scores with the mouse transcriptome for M1 and M2 genes

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_M1_and_M2_concordance_scores.pdf')

par(mar=c(5.1, 4.7, 4.1, 2.1))

# M1
mm = match(genes_upregulated_in_M1_and_downregulated_in_M2$Ensembl, rr$Ensembl)
if (length(which(is.na(mm))) > 0){ genes_upregulated_in_M1_and_downregulated_in_M2 = genes_upregulated_in_M1_and_downregulated_in_M2[-which(is.na(mm)),] }
mm = match(genes_upregulated_in_M1_and_downregulated_in_M2$Ensembl, rr$Ensembl)
rr_M1 = rr[mm,]
disco.score = rr_M1$log2FoldChange * log2(genes_upregulated_in_M1_and_downregulated_in_M2$FC_M1vsM0) * abs(log10(rr_M1$padj) + log10(genes_upregulated_in_M1_and_downregulated_in_M2$p.value_M1vsM0))
oo=order(disco.score)

length(which(disco.score>0))/length(disco.score) # 0.4146341

plot(disco.score[oo],1:length(disco.score[oo]),type='n', axes=FALSE, ylab=NA, xlab='Concordance score', main='Concordance between M1 genes')
axis(1)
par(las=2)
axis(2, at=1:length(disco.score[oo]),labels=rr_M1$gene,cex.axis=0.7)
box(lty='solid')
for(i in 1:length(disco.score[oo])){abline(i,0, col = "lightgray", lty = "dotted")}
points(disco.score[oo],1:length(disco.score[oo]), pch=23, col='black', bg='grey')
abline(v=0, col='black',lty='solid',lwd=0.7)

# M2
mm = match(genes_upregulated_in_M2_and_downregulated_in_M1$Ensembl, rr$Ensembl)
if (length(which(is.na(mm))) > 0){ genes_upregulated_in_M2_and_downregulated_in_M1 = genes_upregulated_in_M2_and_downregulated_in_M1[-which(is.na(mm)),] }
mm = match(genes_upregulated_in_M2_and_downregulated_in_M1$Ensembl, rr$Ensembl)
rr_M2 = rr[mm,]
disco.score = rr_M2$log2FoldChange * log2(genes_upregulated_in_M2_and_downregulated_in_M1$FC_M2vsM0) * abs(log10(rr_M2$padj) + log10(genes_upregulated_in_M2_and_downregulated_in_M1$p.value_M2vsM0))
oo=order(disco.score)

length(which(disco.score>0))/length(disco.score) # 0.3939394

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

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_goseq_pwf_plot.pdf')
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

length(which(GO.KEGG.wall$over_represented_FDR<=0.1)) # 1
GO.KEGG.wall_sig = GO.KEGG.wall[which(GO.KEGG.wall$over_represented_FDR<=0.1),]
# > GO.KEGG.wall_sig
#category over_represented_pvalue under_represented_pvalue numDEInCat
#GO:0010273            1.023194e-06                        1          5
#numInCat                         term ontology over_represented_FDR
#13 detoxification of copper ion       BP           0.01924526

write.csv(GO.KEGG.wall_sig, file = '/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Male_vs_Female_WCM_dataset_DEbySex_ORA.csv')

### make TPM plots of relevant genes

genes_to_plot=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/genes_of_interest.csv', header=TRUE)

# add_gene = data.frame(ensembl=c('ENSG00000154277'), info=c('UCHL1'))
# genes_to_plot=rbind(genes_to_plot, add_gene)

mm=match(genes_to_plot$ensembl, rr$Ensembl)
if(length(which(is.na(mm))) != 0) { genes_to_plot=genes_to_plot[-which(is.na(mm)),] }
mm=match(genes_to_plot$ensembl, rr$Ensembl)
genes_of_interest_rr_subset = rr[mm,]
genes_of_interest_rr_subset$info = genes_to_plot$info
save(genes_of_interest_rr_subset, file = '/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_genes_of_interest_rr_subset.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_genes_of_interest_rr_subset.rda') # genes_of_interest_rr_subset

nrow(genes_to_plot) == length(unique(genes_to_plot$ensembl)) #TRUE

mm = match(genes_to_plot$ensembl, genes$gene_id)
which(is.na(mm)) # none
genes_to_plot$gene_symbol=genes$gene_name[mm]

log2_tpm_plus_one=log2(tpm+1)

all(colnames(log2_tpm_plus_one)==as.vector(ddsClean$names)) #TRUE

mycol=as.vector(ddsClean$Sex)
mycol[which(mycol=="Male")]='blue'
mycol[which(mycol=="Female")]='red'

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_genes_of_interest_TPMs.pdf')

for (i in 1: length(genes_to_plot$ensembl)){

	par(mar=c(5.1,5.3,4.1,2.1))

	index=which(rr$Ensembl==genes_to_plot$ensembl[i])

	zag1=paste0(rr$gene[index],' (', genes_to_plot$info[i], ')')
	zag2=as.expression(bquote(log[2]~"FC" == .(signif(rr$log2FoldChange[index],2))))
	zag3=paste0("FDR = ",signif(rr$padj[index],2))	

	log2_tpm_plus1_subset=as.vector(log2_tpm_plus_one[which(rownames(log2_tpm_plus_one)==genes_to_plot$ensembl[i]),])
	x=factor(as.vector(ddsClean$Sex),levels=c('Male','Female'))

	boxplot(as.vector(log2_tpm_plus1_subset)~x, , xlab='Sex', ylab= as.expression(bquote(log[2]~"(TPM+1)")), main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(as.vector(log2_tpm_plus1_subset) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=c(zag2,zag3), bty='n')

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

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_volcano_plot_DEbySex.pdf')

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
	title='Male vs. female PNET samples',
	subtitle='Discovery dataset',
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
	legendIconSize=15,
	labCol='black',
	labFace='bold'
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

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, ddsClean, rr, vsd

library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

# make TPM matrix
mm=match(rownames(tpm),genes$gene_id)
which(is.na(mm)) # none
rownames(tpm)=genes$gene_name[mm]

# make TPM matrix with removed duplicated genes
tpm_dup_rm=tpm[-c(which(duplicated(rownames(tpm))==TRUE),which(duplicated(rownames(tpm),fromLast=TRUE)==TRUE)),]

## quanTIseq
res_quantiseq = deconvolute(tpm, method="quantiseq", tumor = TRUE)

## EPIC
res_epic = deconvolute(tpm_dup_rm, method="epic", tumor = TRUE)

## CIBERSORT (abs. mode)
set_cibersort_binary("/athena/masonlab/scratch/users/nai2008/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("/athena/masonlab/scratch/users/nai2008/CIBERSORT/LM22.txt")

res_cibersort = deconvolute(tpm_dup_rm, method="cibersort_abs")

## MCP-counter
res_mcpCounter = deconvolute(tpm, method="mcp_counter")

## xCell
res_xCell = deconvolute(tpm, method='xcell')

## make BAR plots comparing cell type fractions

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_inSilico_deconvolution.pdf')

ggplot = function(...) { 	ggplot2::ggplot(...) + theme_bw() }

#quantiseq

res_quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (quanTIseq)')

quantiseq_input = res_quantiseq %>% gather(sample, fraction, -cell_type) 
mm=match(quantiseq_input$sample, ddsClean$names)
quantiseq_input$sex=ddsClean$Sex[mm]

quantiseq_input_male=quantiseq_input[which(quantiseq_input$sex=='Male'),]
quantiseq_input_female=quantiseq_input[which(quantiseq_input$sex=='Female'),]

quantiseq_input_male %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(quantiseq_input_male))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (quanTIseq); Male PNET Samples') +
    theme(legend.position = "none")

quantiseq_input_female %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(quantiseq_input_female))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (quanTIseq); Female PNET Samples') +
    theme(legend.position = "none")

#epic

res_epic %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_epic))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); PNET Samples')

epic_input = res_epic %>% gather(sample, fraction, -cell_type) 
mm=match(epic_input$sample, ddsClean$names)
epic_input$sex=ddsClean$Sex[mm]

epic_input_male=epic_input[which(epic_input$sex=='Male'),]
epic_input_female=epic_input[which(epic_input$sex=='Female'),]

epic_input_male %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(epic_input_male))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); Male PNET Samples') +
    theme(legend.position = "none")

epic_input_female %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(epic_input_female))) +
    ylab('Cell fraction') +
    xlab('Sample') +
    guides(fill=guide_legend(title="Cell Types")) +
    ggtitle('Immune cell type deconvolution (EPIC); Female PNET Samples') +
    theme(legend.position = "none")

dev.off()

## make BOX plots comparing cell type abundance

# quantiseq
res_quantiseq = as.data.frame(res_quantiseq)

cell_type = as.vector(res_quantiseq$cell_type)
abundance_matrix = res_quantiseq[,2:ncol(res_quantiseq)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Sex)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Male')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Female')]),exact=FALSE)
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_deconvolution_cellType_comparisons_quantiseq.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	Sex=factor(colnames(abundance_matrix),levels=c('Male','Female'))

	mycol=as.vector(Sex)
	mycol[which(mycol=="Male")]='blue'
	mycol[which(mycol=="Female")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]

	boxplot(cell_type_data~Sex, , xlab='Sex', ylab= 'Cell type abundance (arbitrary units)', main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(Sex), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

# epic
res_epic = as.data.frame(res_epic)

cell_type = as.vector(res_epic$cell_type)
abundance_matrix = res_epic[,2:ncol(res_epic)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Sex)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Male')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Female')]),exact=FALSE)
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_deconvolution_cellType_comparisons_epic.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	Sex=factor(colnames(abundance_matrix),levels=c('Male','Female'))

	mycol=as.vector(Sex)
	mycol[which(mycol=="Male")]='blue'
	mycol[which(mycol=="Female")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]

	boxplot(cell_type_data~Sex, , xlab='Sex', ylab= 'Cell type abundance (arbitrary units)', main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(Sex), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#CIBERSORT (abs. mode)
res_cibersort = as.data.frame(res_cibersort)

cell_type = as.vector(res_cibersort$cell_type)
abundance_matrix = res_cibersort[,2:ncol(res_cibersort)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Sex)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Male')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Female')]),exact=FALSE)
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_deconvolution_cellType_comparisons_cibersort.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	Sex=factor(colnames(abundance_matrix),levels=c('Male','Female'))

	mycol=as.vector(Sex)
	mycol[which(mycol=="Male")]='blue'
	mycol[which(mycol=="Female")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]

	boxplot(cell_type_data~Sex, , xlab='Sex', ylab= 'Cell type abundance (arbitrary units)', main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(Sex), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#MCP-counter
res_mcpCounter = as.data.frame(res_mcpCounter)

cell_type = as.vector(res_mcpCounter$cell_type)
abundance_matrix = res_mcpCounter[,2:ncol(res_mcpCounter)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Sex)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Male')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Female')]),exact=FALSE)
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_deconvolution_cellType_comparisons_mcpCounter.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	Sex=factor(colnames(abundance_matrix),levels=c('Male','Female'))

	mycol=as.vector(Sex)
	mycol[which(mycol=="Male")]='blue'
	mycol[which(mycol=="Female")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]

	boxplot(cell_type_data~Sex, , xlab='Sex', ylab= 'Cell type abundance (arbitrary units)', main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(Sex), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()

#xCell
res_xCell = as.data.frame(res_xCell)

cell_type = as.vector(res_xCell$cell_type)
abundance_matrix = res_xCell[,2:ncol(res_xCell)]
all(ddsClean$names == colnames(abundance_matrix)) # TRUE
colnames(abundance_matrix) = as.vector(ddsClean$Sex)

p.wilcox.test=vector()
for (i in 1:nrow(abundance_matrix)){
	p=wilcox.test(as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Male')]), as.numeric(abundance_matrix[i,which(colnames(abundance_matrix) == 'Female')]),exact=FALSE)
	p.wilcox.test[i]=p$p.value
}

FDR=p.adjust(p.wilcox.test,method='BH')

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_WCM_dataset_deconvolution_cellType_comparisons_xCell.pdf')

for (i in 1:nrow(abundance_matrix)){

	cell_type_data=as.numeric(abundance_matrix[i,])

	Sex=factor(colnames(abundance_matrix),levels=c('Male','Female'))

	mycol=as.vector(Sex)
	mycol[which(mycol=="Male")]='blue'
	mycol[which(mycol=="Female")]='red'

	par(mar=c(5.1,5.3,4.1,2.1))

	zag1=cell_type[i]

	boxplot(cell_type_data~Sex, , xlab='Sex', ylab= 'Cell type abundance (arbitrary units)', main=zag1, cex.main=1, cex.lab=1.5, cex.axis=1, outline=FALSE, col='lightgrey')
	points(cell_type_data ~ jitter(as.numeric(Sex), amount=0.2), pch =21, col='black', bg=mycol, cex=1)
	legend(x='topright',legend=paste0('FDR = ',signif(FDR[i],2)), bty='n')

}

dev.off()





