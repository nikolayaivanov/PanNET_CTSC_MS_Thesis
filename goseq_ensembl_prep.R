###################################################
## compute median transcript length for each gene
###################################################

library(AnnotationHub)
hub = AnnotationHub()
dm = query(hub, c("EnsDb", "sapiens", "97")) # # use Ensembl - Homo sapiens - release 97 for annotation, to be consistent
edb = dm[["AH73881"]]
tx = transcripts(edb)
tx_length = lengthOf(edb, of = "tx") # details: https://rdrr.io/bioc/ensembldb/man/EnsDb-lengths.html

# convert transcript name to corresponding gene
gene_tx_map = transcripts(edb, columns = c("tx_id", "gene_id"), return.type="data.frame")
mm=match(names(tx_length), gene_tx_map$tx_id)
names(tx_length) = gene_tx_map$gene_id[mm]

tx_lengths_split = split(tx_length, names(tx_length))
median_tx_lengths = unlist(lapply(tx_lengths_split,median))
median_tx_lengths = data.frame(gene_EnsemblID = names(median_tx_lengths), median_length = as.numeric(median_tx_lengths)) # use this to construct the bias.data argumet in the nullp goseq function

###################################################
## make a list where the names are genes and each entry is a vector containing GO categories associated with that gene
###################################################

library(biomaRt)
mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

biomart_GO_results = getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"), mart = mart)

# remove genes with no associated GO term
biomart_GO_results = biomart_GO_results[-which(biomart_GO_results$go_id==""),]

genes_GO_list=split(biomart_GO_results[,1:2],biomart_GO_results$ensembl_gene_id)
gene2cat_GO = lapply(genes_GO_list, '[',,2) # can use this as the gene2cat argument in the goseq function (for GO over-representation analysis)

GO_genes_list=split(biomart_GO_results[,1:2], biomart_GO_results$go_id)
cat2gene_GO=lapply(GO_genes_list, '[',,1)

###################################################
## make a list where the names are genes and each entry is a vector containing KEGG categories associated with that gene
###################################################

biomart_KEGG_results = getBM(attributes = c("ensembl_gene_id", "kegg_enzyme"), mart = mart)

# remove genes with no associated KEGG term
biomart_KEGG_results = biomart_KEGG_results[-which(biomart_KEGG_results$kegg_enzyme==""),]

# separate KEGG enzyme from KEGG pathway IDs
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)
biomart_KEGG_results$kegg_enzyme = ss(biomart_KEGG_results$kegg_enzyme, '\\+', 1)
biomart_KEGG_results$kegg_enzyme = paste0('path:hsa', biomart_KEGG_results$kegg_enzyme)

genes_KEGG_list=split(biomart_KEGG_results,biomart_KEGG_results$ensembl_gene_id)
gene2cat_KEGG = lapply(genes_KEGG_list, '[',,2) # can use this as the gene2cat argument in the goseq function (for KEGG pathway over-representation analysis)

KEGG_genes_list=split(biomart_KEGG_results,biomart_KEGG_results$kegg_enzyme)
cat2gene_KEGG=lapply(KEGG_genes_list, '[',,1)

###################################################
## make a master gene2cat list w/ both GO and KEGG terms
###################################################

mm=match(names(gene2cat_KEGG), names(gene2cat_GO))
length(which(is.na(mm))) # 0

length(names(gene2cat_GO)) == length(unique(names(gene2cat_GO))) #TRUE
length(names(gene2cat_KEGG)) == length(unique(names(gene2cat_KEGG))) #TRUE

for (i in 1:length(gene2cat_GO)){

	mm=match(names(gene2cat_GO[i]),names(gene2cat_KEGG))

	if(!is.na(mm)){

		if (length(mm) == 1){
			if(names(gene2cat_GO[i]) == names(gene2cat_KEGG[mm])){ 
				gene2cat_GO[[i]]=c(gene2cat_GO[[i]],gene2cat_KEGG[[mm]])
			} else if (names(gene2cat_GO[i]) != names(gene2cat_KEGG[mm])){
				print('error: gene names do not match')
			}
		} else if (length(mm) > 1){
			print('error: mm > 1')
		}

	}
}

gene2cat_GOandKEGG = gene2cat_GO

### 

KEGG_term_names = read.table("http://rest.kegg.jp/list/pathway/hsa", sep = "\t", quote = "\"", fill = TRUE, comment.char = "")
colnames(KEGG_term_names) = c('KEGG_ID', 'KEGG_term')

save(median_tx_lengths, gene2cat_GOandKEGG, KEGG_term_names, cat2gene_GO, cat2gene_KEGG, file = '/athena/masonlab/scratch/users/nai2008/items_for_goseq_analysis.rda')
# load('/athena/masonlab/scratch/users/nai2008/items_for_goseq_analysis.rda') # gene2cat_GOandKEGG, KEGG_term_names, median_tx_lengths, cat2gene_GO, cat2gene_KEGG




















