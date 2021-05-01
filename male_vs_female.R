#
phenoData_WCM=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_WCM_PNET_RNAseq_data/PNET_metadata.csv',header=TRUE)
table(phenoData_WCM$Sex)
 # F  M
 # 9 12

phenoData_Chan=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Chan_et_al_GSE118014_PNET_RNAseq_data/pheno_data_for_publication.csv',header=TRUE)
table(phenoData_Chan$Sex)
# Female   Male
#      7      7

phenoData_Fadista=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/_Fadista_GSE50244_control_RNAseq_data/pheno_data_for_publication.csv',header=TRUE)
table(phenoData_Fadista$Sex)
# Female   Male
#     35     54

TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DESeq2_DEbySex_BUNDLE.rda')#  # se, gse, tpm, ddsClean, rr, vsd
rr=rr[-which(is.na(rr$padj)),]
nrow(rr) # 31579
length(which(rr$padj<=0.1)) # 559
m=match(rr$Ensembl,TFs$Ensembl_ID)
length(which(!is.na(m))) # 2734

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Chan_dataset_DESeq2_DEbySex_BUNDLE.rda')
rr=rr[-which(is.na(rr$padj)),]
nrow(rr) # 29574
length(which(rr$padj<=0.1)) # 299
m=match(rr$Ensembl,TFs$Ensembl_ID)
length(which(!is.na(m))) # 2711

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DESeq2_DEbySex_BUNDLE.rda') # se, gse, tpm, ddsClean, rr, vsd
rr=rr[-which(is.na(rr$padj)),]
nrow(rr) # 28753
length(which(rr$padj<=0.1)) # 173
m=match(rr$Ensembl,TFs$Ensembl_ID)
length(which(!is.na(m))) # 2694

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_WCM_dataset_DE_genes_bySex.rda')
DE_genes_WCM=out
nrow(DE_genes_WCM) # 559
length(which(DE_genes_WCM$TF==TRUE)) #55
length(DE_genes_WCM$Ensembl) == length(unique(DE_genes_WCM$Ensembl)) #TRUE

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Chan_dataset_DE_genes_bySex.rda')
DE_genes_Chan=out
nrow(DE_genes_Chan) # 299
length(which(DE_genes_Chan$TF==TRUE)) # 20
length(DE_genes_Chan$Ensembl) == length(unique(DE_genes_Chan$Ensembl)) #TRUE

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Male_vs_Female_Fadista_dataset_DE_genes_bySex.rda')
DE_genes_Fadista=out
nrow(DE_genes_Fadista) # 173
length(which(DE_genes_Fadista$TF==TRUE)) # 15
length(DE_genes_Fadista$Ensembl) == length(unique(DE_genes_Fadista$Ensembl)) #TRUE

mm=match(DE_genes_WCM$Ensembl,DE_genes_Chan$Ensembl)
length(which(!is.na(mm))) # 42

mm=match(DE_genes_Chan$Ensembl,DE_genes_Fadista$Ensembl)
length(which(!is.na(mm))) # 27

mm=match(DE_genes_WCM$Ensembl,DE_genes_Fadista$Ensembl)
length(which(!is.na(mm))) # 30

length(Reduce(intersect, list(DE_genes_WCM$Ensembl,DE_genes_Chan$Ensembl,DE_genes_Fadista$Ensembl))) # 18

# make triple Venn Diagram
library(VennDiagram)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/pdfs/Male_vs_Female_Chan_WCM_Fadista_datasets_TripleVennDiagram.pdf')

overrideTriple=TRUE

draw.triple.venn(area1=559, area2=229, area3=173, n12=42, n23=27, n13=30, n123=19, category =c('PNET\nDiscovery\nDataset','PNET\nValidation\nDataset','Control Dataset'),
    rotation = 1, reverse = FALSE, euler.d =
    TRUE, scaled = TRUE, lwd = rep(2, 3), lty =
    rep("solid", 3), col = rep("black", 3), fill = c('red','blue','green'),
    alpha = rep(0.5, 3), label.col = rep("black", 7), cex
    = c(2,2,2,2,1,2,2), fontface = rep("plain", 7), fontfamily =
    rep("sans", 7), cat.pos = c(-40, 40, 180), cat.dist =
    c(0.07, 0.07, 0.025), cat.col = rep("black", 3),
    cat.cex = rep(.7, 3), cat.fontface = rep("plain", 3),
    cat.fontfamily = rep("sans", 3), cat.just =
    list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
    = "outer", cat.prompts = FALSE, rotation.degree = 0,
    rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
    0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
    sigdigs = 3, direct.area = FALSE, area.vector = 0, margin =0.02)

grid.newpage()

overrideTriple=TRUE

draw.triple.venn(area1=559, area2=229, area3=173, n12=42, n23=27, n13=30, n123=19, category = NA,
    rotation = 1, reverse = FALSE, euler.d = TRUE, scaled = TRUE, lwd = rep(2, 3), lty =
    rep("solid", 3), col = rep("black", 3), fill = c('red','blue','green'),
    alpha = rep(0.5, 3), label.col = rep("black", 7), cex
    = c(2,2,2,2,1,2,2), fontface = rep("plain", 7), fontfamily =
    rep("sans", 7), cat.pos = c(-40, 40, 180), cat.dist =
    c(0.07, 0.07, 0.025), cat.col = rep("white", 3),
    cat.cex = rep(1, 3), cat.fontface = rep("bold", 3),
    cat.fontfamily = rep("sans", 3), cat.just =
    list(c(0.5, 1), c(0.5, 1), c(0.5, 0)), cat.default.pos
    = "outer", cat.prompts = FALSE, rotation.degree = 0,
    rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist =
    0.05, offset = 0, cex.prop = NULL, print.mode = "raw",
    sigdigs = 3, direct.area = FALSE, area.vector = 0, margin =0.02)

dev.off()

# Make a table of genes DE by sex in both the validation and discovery PNET cohorts, and NOT DE by sex in controls
mm=match(DE_genes_WCM$Ensembl,DE_genes_Chan$Ensembl)
validated_genes=DE_genes_WCM[which(!is.na(mm)),1:3]

mm=match(validated_genes$Ensembl,DE_genes_Fadista$Ensembl)
validated_genes_not_DEbysex_in_controls = validated_genes[which(is.na(mm)),]

m=match(validated_genes_not_DEbysex_in_controls$Ensembl,DE_genes_WCM$Ensembl)
validated_genes_not_DEbysex_in_controls$log2FoldChange_discovery = DE_genes_WCM$log2FoldChange[m]
validated_genes_not_DEbysex_in_controls$FDR_discovery = DE_genes_WCM$FDR[m]

m=match(validated_genes_not_DEbysex_in_controls$Ensembl,DE_genes_Chan$Ensembl)
validated_genes_not_DEbysex_in_controls$log2FoldChange_validation = DE_genes_WCM$log2FoldChange[m]
validated_genes_not_DEbysex_in_controls$FDR_validation = DE_genes_Chan$FDR[m]

validated_genes_not_DEbysex_in_controls=validated_genes_not_DEbysex_in_controls[order(validated_genes_not_DEbysex_in_controls$log2FoldChange_discovery),]

validated_genes_not_DEbysex_in_controls_concordantLFC_only=validated_genes_not_DEbysex_in_controls[which(sign(validated_genes_not_DEbysex_in_controls$log2FoldChange_discovery)==sign(validated_genes_not_DEbysex_in_controls$log2FoldChange_validation)),]

write.csv(validated_genes_not_DEbysex_in_controls_concordantLFC_only, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Male_vs_Female_validated_genes_DEbySex_in_PNETs_but_not_in_controls_concordantLFC_only.csv', row.names=FALSE)

# NAI








