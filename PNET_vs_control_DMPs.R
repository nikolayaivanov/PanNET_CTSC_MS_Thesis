##
library(minfi)
library(limma)
library(sva)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/GRset_functionalNorm_SNPs_and_XY_removed.rda.rda')

# library(AnnotationHub)
# hub = AnnotationHub()
# dm = query(hub, c("gencode", "sapiens", "v35"))
# genes=dm[['AH75183']] # AH75183 | Annotated genes for Gencode v31 on hg19 coordinates

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda') # txdb.gencode.v35.hg19, genomicState.gencode.v35.hg19, genes.gencode.v35.hg19, genes.df.gencode.v35.hg19
genes=genes.gencode.v35.hg19

# find DMPs
all_DMPs_rda_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_all_DMPs.rda"
DMPs_proximal_to_genes_rda_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_DMPs_proximal_to_genes.rda"
DMPs_proximal_to_genes_csv_filename="/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/data_tables/PNETvsControls_Chan_and_Syed_DMPs_proximal_to_genes.csv"
pdf_filename='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/pdfs/PNETvsControls_Chan_and_Syed_DMPs.pdf'

pd=pData(GRset)
beta=getBeta(GRset)

dx=factor(as.vector(GRset$Dx), levels=c('Control','PNET')) # Control=0; PNET=1

mod = model.matrix(~dx)

svaobj = sva(beta, mod)

save(svaobj, file='/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_svaobj.rda')
#load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_svaobj.rda')

mod=cbind(mod,svaobj$sv)

probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

probeFit=data.frame(probe=rownames(eb$p.value),
	intercept=probe_fit$coefficients[,1],
	slope=probe_fit$coefficients[,2], 
	p.value=eb$p.value[,2],
	fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr'),
	t=eb$t[,2])

rownames(probeFit)=NULL

o=order(probeFit$fdr)
probeFit=probeFit[o,]
beta=beta[o,]

num_dmps=length(which(probeFit$fdr<=0.05)) # 356,952

percent_dmps=length(which(probeFit$fdr<=0.05))/nrow(probeFit) # 0.7811997

dmps=probeFit[which(probeFit$fdr<=0.05),]

# Determine nearest genes to the DMPs

anno=getAnnotation(GRset)
mm=match(rownames(beta),rownames(anno))
anno=anno[mm,]

rti=as.data.frame(table(factor(anno$Relation_to_Island)))
num_shore_probes=sum(rti$Freq[grep('Shore', rti$Var1, ignore.case=TRUE)])
num_shelf_probes=sum(rti$Freq[grep('Shelf', rti$Var1, ignore.case=TRUE)])
num_island_probes=rti$Freq[grep('Island', rti$Var1, ignore.case=TRUE)]
num_opensea_probes=rti$Freq[grep('OpenSea', rti$Var1, ignore.case=TRUE)]
num_promoter_probes=length(grep('Promoter',anno$Regulatory_Feature_Group, , ignore.case=TRUE))
num_enhancer_probes=length(which(anno$Enhancer=="TRUE"))

anno_DMPs=anno[1:nrow(dmps),]

ad_df=as.data.frame(anno_DMPs)

gr=GRanges(seqnames=ad_df$chr, 
	ranges=IRanges(ad_df$pos, ad_df$pos),
	strand=ad_df$strand)

nn=nearest(gr, genes, select='arbitrary', ignore.strand=TRUE)
hits=genes[nn,]
dist=distance(gr,hits,ignore.strand=TRUE)

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(dmps))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=as.data.frame(cbind(dmps, anno_DMPs[,1:2], distance_to_nearest_gene=dist, Enhancer=Enhancer, Promoter=Promoter))
out$distance_to_nearest_gene=as.numeric(out$distance_to_nearest_gene)

all_dmps=as.data.frame(out)

save(all_dmps, file=all_DMPs_rda_filename) # all DMPs

out1=out[which(out$distance_to_nearest_gene <= 5000 | out$Enhancer == 'TRUE' | out$Promoter == 'TRUE'),]

# Find genes that overlap or are <= 5kb from genes
oo=as.data.frame(findOverlaps(gr, genes, maxgap=5000, select='all', ignore.strand=TRUE))
hits=genes[oo$subjectHits,]
dist=distance(gr[oo$queryHits],hits,ignore.strand=TRUE)

pt1=dmps[oo$queryHits,]
pt2=anno_DMPs[,1:2]
pt2=pt2[oo$queryHits,]

out2=cbind(pt1, pt2, gene=as.vector(hits$Gene), Ensembl=as.vector(hits$Geneid), distance=dist)
save(out2, file=DMPs_proximal_to_genes_rda_filename)

write.csv(out2,file=DMPs_proximal_to_genes_csv_filename, row.names=FALSE) # only the DMPs that overlap or are within 5kb of genes

return=c(num_dmps=nrow(dmps),
	dmps_that_overlap_promoters_or_enhancers=length(which(out$Enhancer=='TRUE' | out$Promoter=='TRUE')),
	dmps_that_overlap_promoters_or_enhancers_or_genes=nrow(out1),
	range_of_abs_values_of_DNAm_differences_in_DMPs_MIN=signif(min(abs(out$slope)),2),
	range_of_abs_values_of_DNAm_differences_in_DMPs_MAX=signif(max(abs(out$slope)),2) )

#                                            num_dmps
#                                         356952.0000
#            dmps_that_overlap_promoters_or_enhancers
#                                         138995.0000
#   dmps_that_overlap_promoters_or_enhancers_or_genes
#                                         341329.0000
# range_of_abs_values_of_DNAm_differences_in_DMPs_MIN
#                                              0.0013
# range_of_abs_values_of_DNAm_differences_in_DMPs_MAX
#                                              0.7300

# output DMP pdf

pdf(file=pdf_filename)

center_hist='TRUE'

if (center_hist=='TRUE'){
	hist(dmps$slope, xlim = c(-max(abs(dmps$slope)),max(abs(dmps$slope))), xlab='PNET versus control\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
		} else {
	hist(dmps$slope, xlab='PNET versus control\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
}

rti_dmps=as.data.frame(table(factor(anno_DMPs$Relation_to_Island)))
Island = rti_dmps$Freq[which(rti_dmps$Var1=='Island')]/num_island_probes
Shore= sum(rti_dmps$Freq[grep('Shore', rti_dmps$Var1, ignore.case=TRUE)])/num_shore_probes
Shelf= sum(rti_dmps$Freq[grep('Shelf', rti_dmps$Var1, ignore.case=TRUE)])/num_shelf_probes
OpenSea= rti_dmps$Freq[which(rti_dmps$Var1=='OpenSea')]/num_opensea_probes
Promoter=length(grep('Promoter',anno_DMPs$Regulatory_Feature_Group, , ignore.case=TRUE))/num_promoter_probes
Enhancer=length(which(anno_DMPs$Enhancer=="TRUE"))/num_enhancer_probes

normalized_Relation_to_Island=data.frame(Island=Island, Shore=Shore, Shelf=Shelf, OpenSea=OpenSea, Promoter=Promoter, Enhancer=Enhancer)

bp=barplot(height=as.matrix(normalized_Relation_to_Island), xlab='Genomic Position', ylab='Normalized Frequency', main='DMPs', col='lightgrey')
labs=signif(as.vector(normalized_Relation_to_Island),3)
text(bp, 0, labs, cex=1, pos=3)

par(mar=c(5.1,5.3,4.1,2.1))

if (num_dmps<100) { n= num_dmps} else {n =100 }

mycol=as.vector(pd$Dx)
mycol=gsub('Control','blue',mycol)
mycol=gsub('PNET','red',mycol)

x=factor(pd$Dx, levels=c('Control','PNET'))

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probe[i])
	#zag2=paste0('p=',signif(probeFit$p.value[i],3), '; t=',signif(probeFit$t[i],3))
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
	boxplot(as.vector(beta[i,])~x,ylab='DNAm (Beta)', xlab='Disease state', main=c(zag1,zag2), cex.main=1, cex.lab=2, cex.axis=2, outline=FALSE, col='lightgrey')
	#stripchart(as.vector(beta[i,])~x, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg=factor(mycol,levels=c('red','blue'))) #bg='gray'
	points(	as.vector(beta[i,]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.2)
}

dev.off()

# DMPs proximal to M1/M2 marker genes

load('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/RNAseq_analysis/rdas/M1_and_M2_gene_sets.rda') #M1_marker_panel,M2_marker_panel,genes_upregulated_in_M1_and_downregulated_in_M2,genes_upregulated_in_M2_and_downregulated_in_M1
load("/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/DNAm_analysis/rdas/PNETvsControls_Chan_and_Syed_DMPs_proximal_to_genes.rda") #out2

nrow(M1_marker_panel) # 20
nrow(M2_marker_panel) # 17

mm=match(M1_marker_panel$Ensembl, out2$Ensembl)
which(is.na) # none!
length(mm) #20

mm=match(M2_marker_panel$Ensembl, out2$Ensembl)
length(which(is.na(mm))) # 1
length(mm[-which(is.na(mm))]) #16















