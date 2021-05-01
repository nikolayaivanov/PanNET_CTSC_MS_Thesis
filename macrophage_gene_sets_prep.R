#
library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97")) # use Ensembl - Homo sapiens - release 97 for annotation, to be consistent
edb = dm[["AH73881"]]
ens.gene.map = genes(edb, columns = c("gene_id", "gene_name"), return.type="data.frame")

M1_marker_panel=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/M1_marker_panel.csv', header=FALSE)
colnames(M1_marker_panel)='Ensembl'

mm=match(M1_marker_panel$Ensembl,ens.gene.map$gene_id)
which(is.na(mm)) #none

M1_marker_panel$gene_symbol=ens.gene.map$gene_name[mm]

##

M2_marker_panel=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/M2_marker_panel.csv', header=FALSE)
colnames(M2_marker_panel)='Ensembl'

mm=match(M2_marker_panel$Ensembl,ens.gene.map$gene_id)
which(is.na(mm)) #none

M2_marker_panel$gene_symbol=ens.gene.map$gene_name[mm]

##

genes_upregulated_in_M1_and_downregulated_in_M2=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/genes_upregulated_in_M1_and_downregulated_in_M2.txt', header=TRUE, sep=" ")

mm=match(toupper(genes_upregulated_in_M1_and_downregulated_in_M2$Gene),ens.gene.map$gene_name)
genes_upregulated_in_M1_and_downregulated_in_M2=genes_upregulated_in_M1_and_downregulated_in_M2[-which(is.na(mm)),]

mm=match(toupper(genes_upregulated_in_M1_and_downregulated_in_M2$Gene),ens.gene.map$gene_name)

genes_upregulated_in_M1_and_downregulated_in_M2$Ensembl=ens.gene.map$gene_id[mm]

gl=split(genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M1_and_downregulated_in_M2$Gene)

for(i in 1:length(gl)){

	if(nrow(gl[[i]]>1)){
		gl[[i]] = gl[[i]][which(gl[[i]]$FC_M1vsM0==max(gl[[i]]$FC_M1vsM0)),]
	}

}

genes_upregulated_in_M1_and_downregulated_in_M2=do.call(rbind.data.frame, gl)

which(duplicated(genes_upregulated_in_M1_and_downregulated_in_M2$Gene)==TRUE) # none

##

genes_upregulated_in_M2_and_downregulated_in_M1=read.table('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/genes_upregulated_in_M2_and_downregulated_in_M1.txt', header=TRUE, sep=" ")

mm=match(toupper(genes_upregulated_in_M2_and_downregulated_in_M1$Gene),ens.gene.map$gene_name)
which(is.na(mm)) #none

genes_upregulated_in_M2_and_downregulated_in_M1$Ensembl=ens.gene.map$gene_id[mm]

gl=split(genes_upregulated_in_M2_and_downregulated_in_M1,genes_upregulated_in_M2_and_downregulated_in_M1$Gene)

for(i in 1:length(gl)){

	if(nrow(gl[[i]]>1)){
		gl[[i]] = gl[[i]][which(gl[[i]]$FC_M1vsM0==max(gl[[i]]$FC_M1vsM0)),]
	}

}

genes_upregulated_in_M2_and_downregulated_in_M1=do.call(rbind.data.frame, gl)

which(duplicated(genes_upregulated_in_M2_and_downregulated_in_M1$Gene)==TRUE) # none

save(M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1, 
		file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1




