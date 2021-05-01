# Discovery dataset: WCM vs Fadista
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/WCM_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/WCM_and_Fadista_datasets_DE_genes_byDx.rda')
out_discovery=out
nrow(out_discovery) # 17,841

# Validation dataset: Chan vs Fadista
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DE_genes_byDx.rda')
out_validation=out
nrow(out_validation) # 17,878

mm=match(out_discovery$Ensembl, out_validation$Ensembl,)
length(which(!is.na(mm))) # 15,207

library(VennDiagram)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/DE_analysis_PNETvsControl_DiscoveryAndValidationDatasetsComparison_VennDiagram.pdf')

overridePairwise=TRUE

draw.pairwise.venn(area1=17841 , area2=17878, cross.area=15207, category = c('Discovery\nAnalysis','Validation\nAnalysis'),
euler.d = TRUE, scaled = TRUE, inverted = FALSE,
ext.text = TRUE, ext.percent = rep(0.05, 3), lwd =
rep(2, 2), lty = rep("solid", 2), col = rep("black",
2), fill = c('blue','green'), alpha = rep(0.5, 2), label.col =
rep("black", 3), cex = rep(1.5, 3), fontface =
rep("plain", 3), fontfamily = rep("sans", 3), cat.pos
= c(-135, 135), cat.dist = rep(0.025, 2), cat.cex =
rep(1, 2), cat.col = rep("black", 2), cat.fontface =
rep("bold", 2), cat.fontfamily = rep("sans", 2),
cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos
= "outer", cat.prompts = TRUE, ext.pos = rep(0, 2),
ext.dist = rep(0, 2), ext.line.lty = "solid",
ext.length = rep(0.95, 2), ext.line.lwd = 1,
rotation.degree = 0, rotation.centre = c(0.5, 0.5),
ind = TRUE, sep.dist = 0.05, offset = 0, cex.prop =
NULL, print.mode = "raw", sigdigs = 3)

dev.off()

########################

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97"))
edb = dm[["AH73881"]]
genes=as.data.frame(genes(edb))

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(DESeq2)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
rr_validation=rr
rr_validation=rr_validation[-which(is.na(rr_validation$padj)),]

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/WCM_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
rr_discovery=rr
rr_discovery=rr_discovery[-which(is.na(rr_discovery$padj)),]

### volcano plot

genes_to_plot=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/genes_of_interest.csv', header=TRUE)

mm = match(genes_to_plot$ensembl, genes$gene_id)
which(is.na(mm)) # none
genes_to_plot$gene_symbol=genes$gene_name[mm]

genes.to.show=as.vector(genes_to_plot[which(genes_to_plot$info %in% c('M1 marker','M2 marker')),])
genes.to.show$col=NA
genes.to.show$col[which(genes_to_plot$info == 'M1 marker')]='blue'
genes.to.show$col[which(genes_to_plot$info == 'M2 marker')]='black'

library(EnhancedVolcano)

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"FC"~"|"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"FC"~"|"~">"~"2"~"&"~'FDR'<='0.1'))

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/DE_analysis_PNETvsControl_DiscoveryAndValidationDatasetsComparison_VolcanoPlot_M1andM2genes.pdf')

# validation analysis

mm=match(genes.to.show$gene_symbol,rr_validation$gene)
genes.to.show_va=genes.to.show[-which(is.na(mm)),]

mm=match(genes.to.show_va$gene_symbol,rr_validation$gene)
rr_validation_subset=rr_validation[mm,]

EnhancedVolcano(toptable = as.data.frame(rr_validation_subset[,c(2,6)]),
	lab = rr_validation_subset$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M1/M2 genes only: Validation analysis',
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
	selectLab=genes.to.show_va$gene_symbol,
	labCol=genes.to.show_va$col,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	#drawConnectors=TRUE
	)

EnhancedVolcano(toptable = as.data.frame(rr_validation_subset[,c(2,6)]),
	lab = rr_validation_subset$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M1/M2 genes only: Validation analysis',
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
	selectLab=genes.to.show_va$gene_symbol,
	labCol=genes.to.show_va$col,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	drawConnectors=TRUE
	)

# discovery analysis

mm=match(genes.to.show$gene_symbol,rr_discovery$gene)
genes.to.show_da=genes.to.show[-which(is.na(mm)),]

mm=match(genes.to.show_da$gene_symbol,rr_discovery$gene)
rr_discovery_subset=rr_discovery[mm,]

EnhancedVolcano(toptable = as.data.frame(rr_discovery_subset[,c(2,6)]),
	lab = rr_discovery_subset$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M1/M2 genes only: Discovery analysis',
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
	selectLab=genes.to.show_da$gene_symbol,
	labCol=genes.to.show_da$col,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	#drawConnectors=TRUE
	)

EnhancedVolcano(toptable = as.data.frame(rr_discovery_subset[,c(2,6)]),
	lab = rr_discovery_subset$gene,
	x = 'log2FoldChange',
	y = 'padj',
	labSize = 3,
	ylab = bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=1,
	title='PNET vs. non-neoplastic pancreatic islets',
	subtitle='M1/M2 genes only: Discovery analysis',
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
	selectLab=genes.to.show_da$gene_symbol,
	labCol=genes.to.show_da$col,
	labFace='bold',
	# labhjust = 0,
	# labvjust = -1,
	drawConnectors=TRUE
	)

dev.off()

## make an output table

out_table=data.frame(gene=genes.to.show$gene_symbol, Ensembl.Gene.ID=genes.to.show$ensembl, gene.category=genes.to.show$info)
write.csv(out_table,file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/M1_and_M2_gene_panels.csv", row.names=FALSE)

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/WCM_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
rr_discovery=rr

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Chan_and_Fadista_datasets_DESeq2_DEbyDx_BUNDLE.rda') # se, gse, tpm, dds, ddsClean, rr, vsd, rld
rr_validation=rr

mm=match(out_table$Ensembl.Gene.ID, rr_discovery$Ensembl)
which(is.na(mm)) # none
out_table$log2FoldChange_discovery=signif(rr_discovery$log2FoldChange[mm],2)
out_table$FDR_discovery=signif(rr_discovery$padj[mm],2)

mm=match(out_table$Ensembl.Gene.ID, rr_validation$Ensembl)
which(is.na(mm)) # none
out_table$log2FoldChange_validation=signif(rr_validation$log2FoldChange[mm],2)
out_table$FDR_validation=signif(rr_validation$padj[mm],2)

write.csv(out_table,file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/DE_analysis_PNETvsControl_DiscoveryAndValidationDatasets_M1andM2genes.csv", row.names=FALSE)

# How many M1 genes are DE by Dx in both discovery and validation datsets?
length(which(out_table$gene.category=='M1 marker'))
length(which(out_table$gene.category=='M1 marker' & out_table$FDR_discovery <=0.1 & out_table$FDR_validation <=0.1))

# How many M2 genes are DE by Dx in both discovery and validation datsets?
length(which(out_table$gene.category=='M2 marker'))
length(which(out_table$gene.category=='M2 marker' & out_table$FDR_discovery <=0.1 & out_table$FDR_validation <=0.1))









