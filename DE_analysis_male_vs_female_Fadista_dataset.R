# Salmon quantification files:
# /athena/masonlab/scratch/users/nai2008/Fadista_etal_control_pancreatic_islets/Salmon_quants_withGCbiasCorrection/

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
library(tximeta)

# Fadista control pancreatic islet data
salmonFiles_Fadista = list.files('/athena/masonlab/scratch/users/nai2008/PNET/Fadista_etal_control_pancreatic_islets/Salmon_quants_withGCbiasCorrection', "quant.sf", recursive=T, full.names=T)
names(salmonFiles_Fadista) = ss(salmonFiles_Fadista,'/',10)
metadata_Fadista_raw=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Fadista_GSE50244_control_RNAseq_data/SraRunTable.csv')
metadata_Fadista=data.frame(Sample_ID=metadata_Fadista_raw$Run, Dx='Control', Loc_or_Met=NA, Sex=metadata_Fadista_raw$gender, Age=metadata_Fadista_raw$Age)

metadata_Fadista$Sex=replace(metadata_Fadista$Sex,which(metadata_Fadista$Sex=='female'),'Female')
metadata_Fadista$Sex=replace(metadata_Fadista$Sex,which(metadata_Fadista$Sex=='male'),'Male')
metadata_Fadista$Sex=factor(metadata_Fadista$Sex, levels=c('Male','Female'))

mm=match(names(salmonFiles_Fadista), metadata_Fadista$Sample_ID)
metadata_Fadista=metadata_Fadista[mm,]
metadata_Fadista=cbind(as.vector(salmonFiles_Fadista),metadata_Fadista)
colnames(metadata_Fadista) = c('files','names','Dx','Loc_or_Met','Sex','Age')

metadata=metadata_Fadista
all(ss(metadata$file,'/',10)==metadata$names) #TRUE

rownames(metadata) = metadata$names

pheno_data_for_publication=metadata[,-1]
colnames(pheno_data_for_publication)[1]='Sample_ID'
write.csv(pheno_data_for_publication, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Fadista_GSE50244_control_RNAseq_data/pheno_data_for_publication.csv', row.names=FALSE)

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
# Number of significant surrogate variables is:  14

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_Fadista_dataset_SV_plots.pdf')

for (i in 2:ncol(svaobj$sv)){
	plot(svaobj$sv[,1],svaobj$sv[,i], col='black', bg=ifelse(metadata$Sex == 'Male', "blue", ifelse(metadata$Sex == 'Female',"red","black")), pch=21, cex=1.5, xlab='SV1',ylab=paste0('SV',i))
	legend('topleft',legend=c('Male','Female'), pch=15, col=c('blue','red'))
}

dev.off()

colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

colData(dds) = as(cbind(as.data.frame(colData(gse)),svaobj$sv),'DataFrame')

design(dds) = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + SV_10 + SV_11 + SV_12 + SV_13 + SV_14 + Sex

# run DE analysis
dds=DESeq(dds)
# 13 rows did not converge in beta
ddsClean <- dds[which(mcols(dds)$betaConv),]
# extract results
rr=results(ddsClean, alpha=0.1, contrast=c('Sex','Male','Female')) 
# contrast = c( the name of a factor in the design formula, name of the numerator level for the fold change, name of the denominator level for the fold change  )
summary(rr)
# out of 30529 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 104, 0.34%
# LFC < 0 (down)     : 69, 0.23%
# outliers [1]       : 0, 0%
# low counts [2]     : 1776, 5.8%
# (mean count < 0)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_Fadista_dataset_Cooks_distances.pdf')

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
save(se, gse, tpm, ddsClean, rr, vsd, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DESeq2_DEbySex_BUNDLE.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, ddsClean, rr, vsd

################ Downstream analysis

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, ddsClean, rr, vsd
rr=rr[-which(is.na(rr$padj)),]

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

### make table of DE genes

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=sig_results$gene, chr=sig_results$chr, Ensembl=sig_results$Ensembl, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

nrow(out) #173

# How many DE genes are TFs?

out$TF=FALSE

TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)

mm=match(out$Ensembl,TFs$Ensembl_ID)
out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) #15

save(out, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DE_genes_bySex.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DE_genes_bySex.rda')

# print table of DE genes
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Male_vs_Female_Fadista_dataset_DE_genes_bySex.csv", row.names=FALSE)

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

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_Fadista_dataset_goseq_pwf_plot.pdf')
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

length(which(GO.KEGG.wall$over_represented_FDR<=0.1)) #0

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

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_Fadista_dataset_volcano_plot_DEbySex.pdf')

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
	title='Male vs. female control pancreatic islet samples',
	subtitle='Control datset',
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




