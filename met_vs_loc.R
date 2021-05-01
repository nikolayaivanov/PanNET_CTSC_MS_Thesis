#
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Met_vs_Loc_WCM_dataset_DE_genes_byMetStatus.rda')
WCM_data=out
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/Met_vs_Loc_Chan_dataset_DE_genes_byMetStatus.rda')
Chan_data=out

mm=match(WCM_data$Ensembl,Chan_data$Ensembl)
length(which(!is.na(mm))) # 34
nrow(WCM_data) # 943
length(which(!is.na(mm)))/nrow(WCM_data) # 0.03605514

WCM_data_validated_subset=WCM_data[-which(is.na(mm)),]
mm=match(WCM_data_validated_subset$Ensembl,Chan_data$Ensembl)
Chan_data_validated_subset=Chan_data[mm,]

colnames(WCM_data_validated_subset)=c('gene','chr','Ensembl','log2FoldChange_WCM','FDR_WCM','TF')
colnames(Chan_data_validated_subset)=c('gene','chr','Ensembl','log2FoldChange_MSK','FDR_MSK','TF')

combined_validated_data = data.frame(gene=WCM_data_validated_subset$gene, chr=WCM_data_validated_subset$chr, Ensembl=WCM_data_validated_subset$Ensembl, 
	log2FoldChange_WCM=WCM_data_validated_subset$log2FoldChange_WCM, FDR_WCM=WCM_data_validated_subset$FDR_WCM, log2FoldChange_MSK=Chan_data_validated_subset$log2FoldChange_MSK, 
	FDR_MSK=Chan_data_validated_subset$FDR_MSK, TF=WCM_data_validated_subset$TF)

write.csv(combined_validated_data, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/data_tables/Met_vs_Loc_DE_genes_byMetStatus_validated.csv', row.names=FALSE)












