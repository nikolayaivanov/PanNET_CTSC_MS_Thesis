##
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# Enable parallelization
library(doParallel)
registerDoParallel(cores = 3)

# Find blocks
pd=pData(GRset)
beta=getBeta(GRset)

dx=factor(as.vector(GRset$Dx), levels=c('Control','PNET')) # Control=0; PNET=1

mod=model.matrix(~dx)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_svaobj.rda')
mod=cbind(mod,svaobj$sv)

cobj=cpgCollapse(GRset, what="Beta")

blocks=blockFinder(cobj$object, design=mod, coef = 2, what = 'Beta', cluster=NULL, nullMethod='bootstrap',
cutoff = 0.1, pickCutoff = FALSE, smooth = TRUE, smoothFunction = locfitByCluster,
B = 500, verbose = TRUE, bpSpan = 2.5 * 10^5)

dat=data.frame(chr=blocks$tab$chr,start=blocks$tab$start,end=blocks$tab$end,p.value=blocks$tab$p.value, FWER=blocks$tab$fwer)

save(blocks, dat, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_blocks.rda')
# load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_blocks.rda')

################ Downstream analysis

library(minfi)
library(rafalib)
library(RColorBrewer)
library(GenomicState)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_blocks.rda')

sig_blocks=blocks$tab[which(blocks$tab$fwer<=0.1),]
nrow(sig_blocks) # 1483

# how many DMRs are hypomethylated in PNETs relative to controls?
length(which(sig_blocks$value<0)) # 1389

#length of the blocks (in bp):
sum(as.integer(sig_blocks$end)-as.integer(sig_blocks$start)+1) # 327,767,206

# make output table with block info, and all the genes they overlap

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda') # txdb.gencode.v35.hg19, genomicState.gencode.v35.hg19, genes.gencode.v35.hg19, genes.df.gencode.v35.hg19
genes=genes.gencode.v35.hg19

out_table=data.frame(chr=sig_blocks$chr, start=sig_blocks$start, end=sig_blocks$end, L=sig_blocks$L, value=sig_blocks$value, p.value=sig_blocks$p.value, fwer=sig_blocks$fwer, overlapping_genes=NA, Ensembl_GeneID=NA)

gr=GRanges(seqnames=sig_blocks$chr, 
	ranges=IRanges(sig_blocks$start, sig_blocks$end))

oo=findOverlaps(gr, genes, maxgap=0, ignore.strand=TRUE)
oo=as.data.frame(oo)

for (i in 1:nrow(sig_blocks)){

	ol=which(oo$queryHits==i)
	index=oo$subjectHits[ol]

	if (length(index)==0){
		out_table$overlapping_genes[i]=NA
		out_table$Ensembl_GeneID[i]=NA
	} else {
		out_table$overlapping_genes[i]=paste(as.vector(genes$Gene[index,]), collapse="; ")
		out_table$Ensembl_GeneID[i]=paste(as.vector(genes$Geneid[index,]), collapse="; ")
	}

}

write.csv(out_table, file="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_DNAm_blocks.csv", row.names=FALSE)

# blocks that overlap M1 genes
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

indexes=vector()
count=0
for (i in 1:nrow(M1_marker_panel)){
  xx=grep(M1_marker_panel$Ensembl[i],out_table$Ensembl_GeneID)
  if(length(xx) !=0){ 
  	indexes=c(indexes,xx)
  	count=count+1
  }
}

indexes=indexes[order(indexes)]
indexes
# 49   79  140  140  872  908 1004

M1_genes_blocks=sig_blocks[indexes,]
nrow(M1_genes_blocks) # 7

nrow(M1_marker_panel) #20
count # 7

out_table_M1_genes=out_table[indexes,]
write.csv(out_table_M1_genes,file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_M1_blocks.csv')

# blocks that overlap M2 genes
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

indexes=vector()
count=0
for (i in 1:nrow(M2_marker_panel)){
  xx=grep(M2_marker_panel$Ensembl[i],out_table$Ensembl_GeneID)
  if(length(xx) !=0){ indexes=c(indexes,xx) 
  count=count+1
  }
}

indexes=indexes[order(indexes)]
indexes
# 10  17 140 158 395

M2_genes_blocks=sig_blocks[indexes,]
nrow(M2_genes_blocks) # 5

nrow(M2_marker_panel) # 17
count # 5

out_table_M2_genes=out_table[indexes,]
write.csv(out_table_M2_genes,file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_M2_blocks.csv')

# plotting

collapsed_CpG_clusters=cpgCollapse(GRset, what='Beta')
cset=collapsed_CpG_clusters$object
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_blocks.rda')

blockPlot(cset=cset,
blocks450=blocks,
coi=pData(GRset)$Dx,
N=length(which(blocks$tab$fwer==0)), #length(which(blocks$tab$fwer<=0.1)) #length(which(blocks$tab$fwer==0))
blockname = "PNET",
filename='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/PNETvsControls_Chan_and_Syed_blocks.pdf',
scale=10,
showMethPanel = TRUE,
showGenePanel=TRUE,
showDiffPanel=TRUE,
showCancerPanel = TRUE,
bty= "o" )

# How many of our sig blocks overlap with cancer DNAm bocks (Hansen et al, 2011, PMID: 21706001)
hansen_blocks=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/cancer_blocks_hansen.csv')

gr=GRanges(seqnames=sig_blocks$chr, 
	ranges=IRanges(sig_blocks$start, sig_blocks$end))

hansen_blocks=GRanges(seqnames=hansen_blocks$Chromosome, 
	ranges=IRanges(hansen_blocks$Start, hansen_blocks$End))

oo=findOverlaps(gr, hansen_blocks, maxgap=0, ignore.strand=TRUE)
oo=as.data.frame(oo)

mm=match(1:nrow(sig_blocks),unique(oo$queryHits))
length(which(is.na(mm))) # 60
# 1423/1483 (96%) of our DNAm blocks overlap with cancer-specific differentially DNA-methylated regions (cDMRs) previously described by Hansen et al, 2011, PMID: 21706001





