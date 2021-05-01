##
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# Enable parallelization
library(doParallel)
registerDoParallel(cores = 4)

# Find bumps
library(bumphunter)

pd=pData(GRset)
dx=factor(as.vector(GRset$Dx), levels=c('Control','PNET')) # Control=0; PNET=1

#Arguments for bumphunter
mod=model.matrix(~dx)
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_svaobj.rda')
mod=cbind(mod,svaobj$sv)
p=getBeta(GRset)

anno=getAnnotation(GRset)
chr=as.vector(anno$chr)
pos=as.vector(anno$pos)

bumps = bumphunterEngine(p, mod, chr = chr, 
pos = pos, cutoff= 0.1, nullMethod = "bootstrap",
smooth=TRUE, B=500)

dat=data.frame(chr=bumps$tab$chr,start=bumps$tab$start,end=bumps$tab$end,p.value=bumps$tab$p.value, FWER=bumps$tab$fwer)

save(bumps, dat, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_DMRs.rda')

################## Downstream analysis

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_DMRs.rda')

library(derfinder)
library(bumphunter)
library(BSgenome.Hsapiens.UCSC.hg19)
library(minfi)
library(devtools)
library(rafalib)
library(RColorBrewer)
library(GenomicState)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# load imprinted genes
load('/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_genes.rda') # ig

# drop imprinted genes without an Ensembl ID
ig=ig[-which(is.na(ig$Ensembl.ID)),]
ig$Ensembl.ID=as.vector(ig$Ensembl.ID)
nrow(ig) #102
length(unique(ig$Ensembl.ID)) #102

# load TFs
TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)
nrow(TFs) # 2765
length(unique(TFs$Ensembl_ID)) # 2765
which(is.na(TFs$Ensembl_ID)) # 0

DMR_table_filename='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_DMRs.csv'

num_dmrs=length(which(dat$FWER<=0.1))
num_dmrs # 4,318

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')
annotation=getAnnotation(GRset)
annotation=as.data.frame(annotation)

bumps$table=bumps$table[order(bumps$table$fwer),]

sigDMRs=bumps$table[bumps$table$fwer<=0.1,]
nrow(sigDMRs) # 4,318

# how many DMRs are hypomethylated in PNETs relative to controls?
length(which(sigDMRs$value<0)) #2078

#length of the DMRs in genomic coordinates:
gl=sum(as.integer(sigDMRs$end)-as.integer(sigDMRs$start)+1) # 160,1919
gl # 160,1919

# make output table with DMR info, and all the genes they overlap
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda') # txdb.gencode.v35.hg19, genomicState.gencode.v35.hg19, genes.gencode.v35.hg19, genes.df.gencode.v35.hg19
genes=genes.gencode.v35.hg19

out_table=data.frame(chr=sigDMRs$chr, start=sigDMRs$start, end=sigDMRs$end, L=sigDMRs$L, genomic_length=gl, value=sigDMRs$value, p.value=sigDMRs$p.value, fwer=sigDMRs$fwer, proximal_genes=NA, Ensembl_GeneID=NA, distance=NA)
# https://support.bioconductor.org/p/90751/

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

oo=findOverlaps(gr, genes, maxgap=5000, ignore.strand=TRUE)
oo=as.data.frame(oo)

oo$dist=distance(gr[oo$queryHits],genes[oo$subjectHits],ignore.strand=TRUE)

for (i in 1:nrow(sigDMRs)){

	ol=which(oo$queryHits==i)

	if( length(ol) != 0 ){
		index=oo$subjectHits[ol]
		out_table$proximal_genes[i]=paste(as.vector(genes$Gene[index,]), collapse="; ")
		out_table$Ensembl_GeneID[i]=paste(as.vector(genes$Geneid[index,]), collapse="; ")
		out_table$distance[i]=paste(oo$dist[ol], collapse="; ")
	} else {
		out_table$proximal_genes[i]='No overlapping or proximal genes (up to 5kb away)'
	}

}

write.csv(out_table, DMR_table_filename, row.names=FALSE)

info=list()
olg=as.vector(genes$Geneid[oo$subjectHits])

# how many overlapping genes are imprinted genes?

mm_ig=match(olg,ig$Ensembl.ID)

if(length(which(!is.na(mm_ig))) !=0 ){

	if (length(which(is.na(mm_ig))) != 0) { mm_ig=mm_ig[-which(is.na(mm_ig))] }
	
	info[[1]]=as.vector(ig$Gene[mm_ig])

} else { info[[1]] = 'None of the DMRs overlap or are within 5kb of imprinted genes'}

names(info)[1]='Imprinted_genes'

length(info[[1]]) # 88

# how many pverlapping genes are TFs?

mm_tf=match(olg,TFs$Ensembl_ID)

if(length(which(!is.na(mm_tf))) !=0 ){

	if ( length(which(is.na(mm_tf))) != 0 ) { mm_tf=mm_tf[-which(is.na(mm_tf))] }
	
	info[[2]]=as.vector(TFs$HGNC_symbol[mm_tf])

} else { info[[2]] = 'None of the DMRs overlap or are within 5kb of TFs'}

names(info)[2]='Transcription_factors'

length(info[[2]]) #802

# match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

match_DMRs_to_genes=matchGenes(gr, genes)

for (i in 1:nrow(match_DMRs_to_genes)){
	if(is.na(match_DMRs_to_genes$name[i])) { match_DMRs_to_genes$name[i]=as.vector(match_DMRs_to_genes$Geneid[i]) } 
}

## plotting DMRs

chr=as.vector(annotation$chr)
pos=as.vector(annotation$pos)
cluster=clusterMaker(chr,pos,maxGap=500)

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/PNETvsControls_Chan_and_Syed_DMRs.pdf')

dmrPlot(regions=sigDMRs,
	p=getBeta(GRset),
	chr=as.vector(annotation$chr),
	pos=as.vector(annotation$pos),
	cluster=cluster,
	genes=match_DMRs_to_genes,
	coi=pData(GRset)$Dx,
	build="hg19",
	number=length(which(bumps$table$fwer==0)),
	species = "human",
	Jitter = TRUE,
	cols=c("royalblue2", "red"),
	lines = FALSE,
	linesSmooth = TRUE,
	title = TRUE,
	Legend = TRUE,
	colorRamp = FALSE, 
	meanSmooth=TRUE,
	plotCpG = TRUE,
	geneAnno = "gene" )

dev.off()

# DMRs that overlap M1 genes
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

indexes=vector()
for (i in 1:nrow(M1_marker_panel)){
  xx=grep(M1_marker_panel$Ensembl[i],out_table$Ensembl_GeneID)
  if(length(xx) !=0){ indexes=c(indexes,xx) }
}

M1_genes_bumps=sigDMRs[indexes,]
nrow(M1_genes_bumps) # 9

out_table_M1_genes=out_table[indexes,]
write.csv(out_table_M1_genes,file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_M1_DMRs.csv')

# match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=M1_genes_bumps$chr, 
	ranges=IRanges(M1_genes_bumps$start, M1_genes_bumps$end))

match_DMRs_to_genes=matchGenes(gr, genes)

for (i in 1:nrow(match_DMRs_to_genes)){
	if(is.na(match_DMRs_to_genes$name[i])) { match_DMRs_to_genes$name[i]=as.vector(match_DMRs_to_genes$Geneid[i]) } 
}

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/PNETvsControls_Chan_and_Syed_M1_DMRs.pdf')

# source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# chr=as.vector(annotation$chr)
# pos=as.vector(annotation$pos)
# cluster=clusterMaker(chr,pos,maxGap=500)

dmrPlot(regions=M1_genes_bumps,
	p=getBeta(GRset),
	chr=as.vector(annotation$chr),
	pos=as.vector(annotation$pos),
	cluster=cluster,
	genes=match_DMRs_to_genes,
	coi=pData(GRset)$Dx,
	build="hg19",
	number=nrow(M1_genes_bumps),
	species = "human",
	Jitter = TRUE,
	cols=c("royalblue2", "red"),
	lines = FALSE,
	linesSmooth = TRUE,
	title = TRUE,
	Legend = TRUE,
	colorRamp = FALSE, 
	meanSmooth=TRUE,
	plotCpG = TRUE,
	geneAnno = "gene" )

dev.off()

# DMRs that overlap M2 genes
load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1

indexes=vector()
for (i in 1:nrow(M2_marker_panel)){
  xx=grep(M2_marker_panel$Ensembl[i],out_table$Ensembl_GeneID)
  if(length(xx) !=0){ indexes=c(indexes,xx) }
}

M2_genes_bumps=bumps$tab[indexes,]

M2_genes_bumps=sigDMRs[indexes,]
nrow(M2_genes_bumps) # 2

out_table_M2_genes=out_table[indexes,]
write.csv(out_table_M2_genes,file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_M2_DMRs.csv')

# match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=M2_genes_bumps$chr, 
	ranges=IRanges(M2_genes_bumps$start, M2_genes_bumps$end))

match_DMRs_to_genes=matchGenes(gr, genes)

for (i in 1:nrow(match_DMRs_to_genes)){
	if(is.na(match_DMRs_to_genes$name[i])) { match_DMRs_to_genes$name[i]=as.vector(match_DMRs_to_genes$Geneid[i]) } 
}

pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/PNETvsControls_Chan_and_Syed_M2_DMRs.pdf')

# chr=as.vector(annotation$chr)
# pos=as.vector(annotation$pos)
# cluster=clusterMaker(chr,pos,maxGap=500)

dmrPlot(regions=M2_genes_bumps,
	p=getBeta(GRset),
	chr=as.vector(annotation$chr),
	pos=as.vector(annotation$pos),
	cluster=cluster,
	genes=match_DMRs_to_genes,
	coi=pData(GRset)$Dx,
	build="hg19",
	number=nrow(M2_genes_bumps),
	species = "human",
	Jitter = TRUE,
	cols=c("royalblue2", "red"),
	lines = FALSE,
	linesSmooth = TRUE,
	title = TRUE,
	Legend = TRUE,
	colorRamp = FALSE, 
	meanSmooth=TRUE,
	plotCpG = TRUE,
	geneAnno = "gene" )

dev.off()





































