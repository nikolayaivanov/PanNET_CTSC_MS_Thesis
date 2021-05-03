#
#source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# Function to capitalize teh first letter
# from https://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

cleaningP = function(y, mod, P=ncol(mod)) {
  X=mod
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

# DMP validation function
validate_dmps = function(GRset.funnorm,studyID) {

GRset.funnorm_subset=GRset.funnorm[,which(pData(GRset.funnorm)$studyID==studyID)]
beta=getBeta(GRset.funnorm_subset)

pd=pData(GRset.funnorm_subset)

mod=model.matrix(~as.factor(pd$sex))
svaobj = sva(beta, mod)

mod = model.matrix(~as.factor(pd$sex))
mod=cbind(mod,svaobj$sv)
probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

slope=probe_fit$coefficients[,2]
intercept=probe_fit$coefficients[,1]
p.value=eb$p.value[,2]
t=eb$t[,2]
fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr')

probeFit=data.frame(probes=rownames(eb$p.value), slope=slope, intercept=intercept, 
  p.value=p.value, fdr=fdr, t=t)
rownames(probeFit)=NULL

file=paste('/athena/masonlab/scratch/users/nai2008/DNAm/rdas/probeFit_',studyID,'.rda',sep='')
save(probeFit,file=file)

dmps=which(probeFit$fdr<=.05)

return(length(dmps))

}

# Gene Ontology functions

dogo_entrez <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as entrez
  cond=FALSE, ontology = "BP"){
    if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
    golib="org.Rn.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Rn.egREFSEQ2EG
  }
  require(GOstats)
  x=names
  x=x[!is.na(x)]

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

dogo_RefSeq <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as RefSeq
	cond=FALSE, ontology = "BP"){
    if(species=="human"){
		golib="org.Hs.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
		golib="org.Mm.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
		golib="org.Rn.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Rn.egREFSEQ2EG
	}
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA)) # convert inputted RefSeq genes (user set) to entrez 
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA)) # convert inputted RefSeq genes (universe) to entrez  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

dogo_GeneSymbols <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as Gene Symbol
  cond=FALSE, ontology = "BP"){
    if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egSYMBOL2EG
  } else  if (species == "mouse") {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egSYMBOL2EG
  } else if (species == "rat") {
    golib="org.Rn.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Rn.egSYMBOL2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA)) # convert inputted Gene Symbols (user set) to entrez 
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA)) # convert inputted Gene Symbols (universe) to entrez  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

# wrapper for string split and sapply
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

# command that will return % of variance explained by each PC
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100

# function that takes 2 arguments x1 and x2, and matches all elements of x1 to ALL elements of x2, and returns an index vector corresponding to x2

match_all = function (x1, x2){

indexes=list()

  for (i in 1:length(x1)){ 

    aa = which( x2 %in% x1[i]) 

    if(length(aa)>0){indexes[[i]] = aa 
      } else if (length(aa)==0) { indexes[[i]] = NA }
  }

oo=unlist(indexes)

return(oo)

}

# function to generate a summary file a from fastqc files

fastqcSummary <- function(files){ 
#files=paths to unzipped fastqc folders 

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

txt_files=paste0(files,'/fastqc_data.txt')

  for (i in 1:length(txt_files)){
  
    scan=as.vector(sapply(txt_files[i], function(x) scan(x,"",sep="\n")))

    ## Total Sequences
    ts=scan[grep('Total Sequences',scan)]
    Total_Sequences=as.numeric(ss(ts,'Sequences\t',2))

    ## Per base sequence quality 
    Per_base_sequence_quality=read.table(txt_files[i],sep='\t',comment.char="",skip=grep('Per base sequence quality',scan)+1, nrows=(grep('Per tile sequence quality',scan)-grep('Per base sequence quality',scan)-3))
    if (ncol(Per_base_sequence_quality)==7){
      colnames(Per_base_sequence_quality)=c('Base', 'Mean', 'Median', 'Lower_Quartile', 'Upper_Quartile', '10th_Percentile', '90th_Percentile') 
    }

    ## Adapter Content
    Adapter_Content=read.table(txt_files[i],sep='\t',comment.char="",skip=grep('Adapter Content',scan)+1, nrows=(grep('Kmer Content',scan)-grep('Adapter Content',scan)-3)) 
    if(ncol(Adapter_Content)==4){
    colnames(Adapter_Content)=c('Position', 'Illumina_Universal_Adapter', 'Illumina_Small_RNA_Adapter', 'Nextera_Transposase_Sequence')
    } else if (ncol(Adapter_Content)==5){
    colnames(Adapter_Content)=c('Position', 'Illumina_Universal_Adapter', 'Illumina_Small_RNA_Adapter', 'Nextera_Transposase_Sequence', 'SOLID_Small_RNA_Adapter')
    }

    fastqcSummary=list(Total_Sequences,Per_base_sequence_quality,Adapter_Content)
    names(fastqcSummary)=c('total_sequences','Per_base_sequence_quality','Adapter_Content')

    output_dir=files[i]
    o=paste0(output_dir,'/fastqcSummary.rda')
    save(fastqcSummary,file=o)

  } 

  #if (length(files)==1) { return(fastqcSummary) }

} #end of function


## finding DMRs

DMRs_v2 = function(input_GRset, DMRs_filename, B) {

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 4)

# Find bumps
library(bumphunter)

pd=pData(input_GRset)
neg_control_PCs=pd[,13:21]
sex=factor(as.vector(pd$sex), levels=c('Female','Male')) # Female = 0; Male=1;

#Arguments for bumphunter
mod=model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)
p=getBeta(input_GRset)

anno=getAnnotation(input_GRset)
chr=as.vector(anno$chr)
pos=as.vector(anno$pos)

bumps = bumphunterEngine(p, mod, chr = chr, 
pos = pos, cutoff= 0.1, nullMethod = "bootstrap",
smooth=TRUE, B=B)

dat=data.frame(chr=bumps$tab$chr,start=bumps$tab$start,end=bumps$tab$end,p.value=bumps$tab$p.value, FWER=bumps$tab$fwer)

save(bumps, dat, file=DMRs_filename)

num_dmrs=length(which(dat$FWER<=0.1))

return(num_dmrs)

}

## plot DMRs
dmrPlot = function(regions, p, chr, pos, cluster, genes, coi, build="hg19", number=100, species = "human",
  Jitter = FALSE, cols=NULL, lines = FALSE, linesSmooth = TRUE, title = TRUE, Legend = TRUE, colorRamp = FALSE, 
  meanSmooth=TRUE, plotCpG = TRUE, geneAnno = "gene") {
  
  load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda') # txdb.gencode.v35.hg19, genomicState.gencode.v35.hg19, genes.gencode.v35.hg19, genes.df.gencode.v35.hg19
  gs=genomicState.gencode.v35.hg19$fullGenome

  #require(bumphunter)
  #require(derfinder)
  gr = GRanges(regions$chr, IRanges(regions$start, regions$end))
  
  if(build == "hg18") {
    cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Hsapiens.UCSC.hg18")
  }
  
  if(build == "hg19") {
    cpg.cur = read.table("http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Hsapiens.UCSC.hg19")
  }
  
  if(build == "mm9") {
    cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-mm9.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Mmusculus.UCSC.mm9")
  }
    
  ocpgi=data.frame(chr=I(cpg.cur[,1]), 
    start=as.numeric(cpg.cur[,2]), 
    end=as.numeric(cpg.cur[,3]))
  ADD1 = 1500; PAD = 10
  
  gr2 = GRanges(regions$chr, IRanges(regions$start - ADD1, regions$end + ADD1))
  anno = annotateRegions(gr2, gs)$annotationList
  
  # check regions
  if(is.numeric(coi)) {
    groups=coi
    gNames= sort(unique(coi))
  }

  if(is.character(coi) | is.factor(coi)) {
    groups = factor(coi)
    gNames= levels(groups)
  }
  gIndexes=split(1:length(groups),groups)
  
  brewer.n=max(3,min(length(unique(coi)),11))

  
  if(is.null(cols)) {
    mypar(brewer.n=brewer.n)
  } else if(length(cols) == 1 & cols[1] %in% rownames(brewer.pal.info)) {
    mypar(brewer.name=cols, brewer.n=brewer.n)
  } else {
    mypar()
    palette(cols)
  } 
    
  if(colorRamp) {
    pal = colorRampPalette(palette()[3:brewer.n])
    palette(pal(brewer.n))
  }
  
  cat("Plotting.\n")
  for(j in 1:(min(nrow(regions), number))) {
  
    layout(matrix(1:2,ncol=1),heights=c(0.7,0.3))

    # first plot, region
    par(mar=c(0,4.5,0.25,1.1),oma=c(0,0,2,0))
    
    Index=(regions[j,7]-PAD):(regions[j,8]+PAD)
    Index = Index[Index %in% seq(along=cluster)]
    Index=Index[cluster[Index]==regions[j,6]]
    
    # make DMR plots
    if(Jitter) {
      posx = matrix(pos[Index], nc = ncol(p), nr = length(Index),
        byrow=  FALSE)
      posx = t(apply(posx, 1, function(x) jitter(x,amount = 12)))
      
    } else posx = pos[Index]
          
    matplot(posx, p[Index,], ylim = c(0,1),
      ylab = "", xlab = "",xaxt = "n",cex=0.7,
      cex.axis = 1.7, cex.lab = 1.7, pch=21,
      bg = as.numeric(factor(groups)),col="black",
      xlim = range(pos[Index]), yaxt="n")
    axis(2, at = c(0.2, 0.5, 0.8), cex.axis=1.7)

    legend("topright", paste0("fwer = ", signif(regions$fwer[j],3)), bty='n')

    xx=pos[Index]
    for(k in seq(along=gIndexes)){
      if(length(gIndexes[[k]]) == 1) yy=p[Index,gIndexes[[k]]]
      if(length(gIndexes[[k]]) > 1)   yy=rowMeans(p[Index,gIndexes[[k]]])
      if(meanSmooth) { 
        fit1=loess(yy~xx,degree=1,span=300/(36*length(xx)),
          family="symmetric")
        lines(xx,fit1$fitted,col=k,lwd=2)
      } else  lines(xx,yy,col=k,lwd=2)

    }
    
    mtext("Methylation",side=2, line = 2.5,cex=1.8)

    if(Legend) {
      if(length(unique(coi)) < 4) {
        legend("topleft",legend=gNames,col=1:length(gNames),
        lty=1, lwd = 4,cex=1,bty="n")
      } else {
        legend("topleft",legend=gNames,col=1:length(gNames),
          pch=19, pt.cex = 2,cex=1.1, nc = 6,bty="n")
      } 
    }
        
    abline(v=(regions$start[j]-15),lty=1)
    abline(v=(regions$end[j]+15),lty=1)

    if(title) mtext(paste0(genes$name[j],"; ", 
      genes$distance[j],"bp from tss:",genes$description[j]), outer=T,cex=1.3)
      
    # add cpgs
    if(plotCpG) {
      thechr=as.character(regions$chr[j])
      chrName = strsplit(thechr, "r")[[1]][2]
      chrName = paste("Chromosome",chrName)
      
      start = pos[Index[1]]
      end = pos[Index[length(Index)]]
      ocpgi2=ocpgi[ocpgi$chr%in%unique(as.character(thechr)),]
      
      ##PLOT CPG ISLANDS
      if(species=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(species=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      
      subseq<-subseq(seq,start=start,end=end)
      cpgs=start(matchPattern("CG",subseq))+start-1

      if(plotCpG) rug(cpgs,col="black")  #  previously had 'Rug'

      Index1 = which(ocpgi2[,1]==as.character(thechr) &
           ((ocpgi2[,2] > start & ocpgi2[,2]< end) |
            (ocpgi2[,3] > start & ocpgi2[,3]< end)))
      if(length(Index1)>0) sapply(Index1,function(j) rug(unlist(ocpgi2[j,2:3]),  #  previously had 'Rug'
             col="darkgreen",lwd=3,side=1))
    }
    
    # plot 3
    ##PLOT GENES
    par(mar=c(3.5,4.5,0.25,1.1))

    plot(0,0, type="n", xlim=range(xx),ylim=c(-1.5,1.5),yaxt="n",ylab="",
       xlab="",cex.axis = 1.5, cex.lab =1.5)
    a = as.data.frame(anno[[j]])
    Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
    Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
    Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
    axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
    abline(h=0,lty=3)
    for(k in 1:nrow(a)) {
      polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
        Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
    }
    
    if(sum(a$theRegion=="exon") > 0) {
      e = a[a$theRegion=="exon",]
      s2 = Strand[a$theRegion=="exon"]
      g = unlist(e$symbol)
      g[is.na(g)] = ""
      if(length(g) > 0) text(x=e$start + e$width/2,
        y=s2*0.75, g,font=1,pos=s2+2,cex=1.2)
    }
        
    mtext("Genes",side=2, line = 2.5,cex=1.5)
    mtext(chrName,side=1, line = 2,cex=1.4)

    abline(v=(pos[regions[j,7]]-15),lty=1)
    abline(v=(pos[regions[j,8]]+15),lty=1)
    
  }

}

## finding DNAm blocks

blocks_v2 = function(input_GRset, blocks_filename, B) {

library(minfi)

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 3)

# Find blocks
pd=pData(input_GRset)
neg_control_PCs=pd[,13:21]
sex=factor(as.vector(pd$sex), levels=c('Female','Male')) # Female = 0; Male=1;

mod=model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)

cobj=cpgCollapse(input_GRset, what="Beta")

blocks=blockFinder(cobj$object, design=mod, coef = 2, what = 'Beta', cluster=NULL, nullMethod='bootstrap',
cutoff = 0.1, pickCutoff = FALSE, smooth = TRUE, smoothFunction = locfitByCluster,
B = B, verbose = TRUE, bpSpan = 2.5 * 10^5)

dat=data.frame(chr=blocks$tab$chr,start=blocks$tab$start,end=blocks$tab$end,p.value=blocks$tab$p.value, FWER=blocks$tab$fwer)

save(blocks, dat, file=blocks_filename)

num_blocks=length(which(dat$FWER<=0.1))

return(num_blocks)

}

## plot blocks

## cset: result of minfi::cpgCollapse()$object
## blocks450: results of blockFinder
## coi = covariate of interest
## N: the number of blocks to plot, default=10
## blockname = name of block track, default='coi'
## filename = where to save plots
## scale, in kb. default = 100

blockPlot = function(cset,
  blocks450,
  coi,
  N=100,
  blockname = "coi",
  filename,
  scale=100,
  showMethPanel = TRUE,
  showGenePanel=TRUE,
  showDiffPanel=TRUE, 
  showCancerPanel = FALSE,
  bty= "o") {

  panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
  
  #require(GenomicRanges)

  blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
  colIds=match(c("chr","start","end"),names(blocks450$table))
  mcols(blocksTable)=blocks450$table[-colIds]

  plotRegion = blocksTable[1:N]

  ## annotation based on ensembl
  cat("Loading Annotation.\n")
  #load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/genomicState.rda")
  load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda') # txdb.gencode.v35.hg19, genomicState.gencode.v35.hg19, genes.gencode.v35.hg19, genes.df.gencode.v35.hg19
  genomicState = genomicState.gencode.v35.hg19
  gs = genomicState$fullGenome
  oo = findOverlaps(blocksTable, gs)
  anno = split(gs[subjectHits(oo)], queryHits(oo))

  # cancer blocks
  if(showCancerPanel) {
    # read in cancer DNAm bocks (Hansen et al, 2011, PMID: 21706001)
    #load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
    blocks=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/cancer_blocks_hansen.csv')
    blocks=GRanges(blocks)
    genome(blocks)="hg19"
    cancerBlocks = blocks
  }

  cat("Ploting.")
  pdf(filename,height=5,width=10)
  par(bty=bty)

  for(i in seq(along=plotRegion)) { #//
    cat(".")
    r = plotRegion[i]
    tmp=subsetByOverlaps(cset,r)
    tmp450=sort(subsetByOverlaps(blocksTable,r))
    if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
                      
    beta=getBeta(tmp)
    x=start(tmp)

    ii=cset %over% r
    d=blocks450$coef[ii]
    sd=blocks450$fitted[ii]
    ## which rows
    Index=split(seq_along(coi),coi)
    mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
    smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

    ## paneling, from hector
    mypar(1,1, brewer.name = "Set1")
    par(mar=par()$mar+c(0,3,0,0))
    omar=par()$mar
    cmar=omar
    cmar[1]=.5
    par(mar=cmar)
    
    mycols=c("royalblue2", "red") ###
    palette(mycols) ###

    layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
    if(showMethPanel) {
      matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
      xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
      matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
      legend("bottomright",col=seq(along=levels(factor(coi))),
        lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
        
      segments(min(x), .2, min(x)+scale*1000, .2)
      text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
      axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
      mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
      mtext("Methylation",side=2, line = 2.5,cex=1.5)

      legend("topright", paste0("fwer = ",
        signif(blocks450$tab$fwer[i],3)),bty=bty)
    }   
    # annotation
    if(showGenePanel) {
      cmar=omar
      if(!is.na(showDiffPanel)) {
        cmar[3]=0.5
      } else {
        cmar[c(1,3)]=c(0.5,0.5)
      }
      par(mar=cmar)

      plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
         xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
      a = as.data.frame(anno[[i]])
      Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
      Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
      Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
      axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
      abline(h=0,lty=3)
      for(k in 1:nrow(a)) {
        polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
          Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
      }

      ## by gene
      g = split(a, sapply(a$symbol,"[", 1))
      # g = split(a, a$Symbol)
      s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
      g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
      
      if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
          
      mtext("Genes",side=2, line = 2.5,cex=1.5)
      if(!showDiffPanel) {
        xtick=pretty(x)
        axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
      }
    }

    ## mean diff
    if(showDiffPanel) {
      cmar=omar
      cmar[c(1,3)]=c(.5,.5)
      par(mar=cmar)

      zz=granges(tmp)

      matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
      axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
      xtick=pretty(x)
      axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
      mtext("Diff",side=2, line = 2.5,cex=1.5)

      ii=which(zz$type=="OpenSea")
      blockgroup=zz$blockgroup[ii]

      blockIndexes=split(seq(along=blockgroup),blockgroup)
      for (ind in blockIndexes) {
        ind=ii[ind]
        lines(x[ind], sd[ind], lwd=2.5,col="black")
      }


      points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
      axis(side=2,at=c(-2,0,2))
      abline(h=0,lty=2,col="black")

      cmar=omar
      cmar[3]=.5
      par(mar=cmar)
      matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
        ylab="",ylim=c(0,2),bty="n")
    }
    
    #  browser()
    if(showCancerPanel) {
      col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
      rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
      if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$Direction.of.Methylation.Change=="hypo","blue","red"))
      axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
      legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
    }
  } #//
  
  dev.off()

}

# function to plot DNAm age (Horvath) wrt actual (chronological) age

DNAmAge_plots = function(females_DNAmAge , males_DNAmAge, pdf_filename) {

pdf(pdf_filename)

## all ages

fit=lm(males_DNAmAge$DNAmAge~males_DNAmAge$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(males_DNAmAge$RealAge, males_DNAmAge$DNAmAge, method='spearman', alternative='two.sided', exact=FALSE)
male_legend_label=paste0('Male: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

par(mar=c(5.1,5,4.1,2.1))

xlim=c(min(c(males_DNAmAge$RealAge,females_DNAmAge$RealAge)),max(c(males_DNAmAge$RealAge,females_DNAmAge$RealAge)))
ylim=c(min(c(males_DNAmAge$DNAmAge,females_DNAmAge$DNAmAge)),max(c(males_DNAmAge$DNAmAge,females_DNAmAge$DNAmAge)))

plot(males_DNAmAge$RealAge, males_DNAmAge$DNAmAge, xlab='Chronological Age',ylab='DNAm Age', xlim=xlim, ylim=ylim,
  pch=21, col='black', bg='blue', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
abline(intercept,slope,lty='solid',col='blue',lwd=4)

fit=lm(females_DNAmAge$DNAmAge~females_DNAmAge$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(females_DNAmAge$RealAge, females_DNAmAge$DNAmAge, method='spearman', alternative='two.sided', exact=FALSE)
female_legend_label=paste0('Female: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

points(females_DNAmAge$RealAge, females_DNAmAge$DNAmAge,pch=21, col='black', bg='red', cex=2)
abline(intercept,slope,lty='solid',col='red',lwd=4)

abline(0,1,lty=2,col='black',lwd=4)

legend('topleft', c(male_legend_label, female_legend_label), cex=.8)
legend('bottomright', as.expression(bquote(bold("All ages"))), cex=1.5, bty='n')
#legend("topleft", c('Males', 'Females','Identity'), lwd=4, col=c('blue','red','black'), lty=c(1,1,2), cex=1.3)

## peds only

females_DNAmAge_peds=females_DNAmAge[which(females_DNAmAge$RealAge<=21),]
males_DNAmAge_peds=males_DNAmAge[which(males_DNAmAge$RealAge<=21),]

fit=lm(males_DNAmAge_peds$DNAmAge~males_DNAmAge_peds$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(males_DNAmAge_peds$RealAge, males_DNAmAge_peds$DNAmAge, method='spearman', alternative='two.sided', exact=FALSE)
male_legend_label=paste0('Male: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

par(mar=c(5.1,5,4.1,2.1))

xlim=c(min(c(males_DNAmAge_peds$RealAge,females_DNAmAge_peds$RealAge)),max(c(males_DNAmAge_peds$RealAge,females_DNAmAge_peds$RealAge)))
ylim=c(min(c(males_DNAmAge_peds$DNAmAge,females_DNAmAge_peds$DNAmAge)),max(c(males_DNAmAge_peds$DNAmAge,females_DNAmAge_peds$DNAmAge)))

plot(males_DNAmAge_peds$RealAge, males_DNAmAge_peds$DNAmAge, xlab='Chronological Age',ylab='DNAm Age', xlim=xlim, ylim=ylim,
  pch=21, col='black', bg='blue', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
abline(intercept,slope,lty='solid',col='blue',lwd=4)

fit=lm(females_DNAmAge_peds$DNAmAge~females_DNAmAge_peds$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(females_DNAmAge_peds$RealAge, females_DNAmAge_peds$DNAmAge, method='spearman', alternative='two.sided', exact=FALSE)
female_legend_label=paste0('Female: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

points(females_DNAmAge_peds$RealAge, females_DNAmAge_peds$DNAmAge,pch=21, col='black', bg='red', cex=2)
abline(intercept,slope,lty='solid',col='red',lwd=4)

abline(0,1,lty=2,col='black',lwd=4)

legend('topleft', c(male_legend_label, female_legend_label), cex=.8)
legend('bottomright', as.expression(bquote(bold("0-21 y.o."))), cex=1.5, bty='n')
#legend("topleft", c('Males', 'Females','Identity'), lwd=4, col=c('blue','red','black'), lty=c(1,1,2), cex=1.3)

## adults only

females_DNAmAge_adults=females_DNAmAge[which(females_DNAmAge$RealAge>21),]
males_DNAmAge_adults=males_DNAmAge[which(males_DNAmAge$RealAge>21),]

fit=lm(males_DNAmAge_adults$DNAmAge~males_DNAmAge_adults$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(males_DNAmAge_adults$RealAge, males_DNAmAge_adults$DNAmAge, method='spearman', alternative='two.sided', exact=FALSE)
male_legend_label=paste0('Male: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

par(mar=c(5.1,5,4.1,2.1))

xlim=c(min(c(males_DNAmAge_adults$RealAge,females_DNAmAge_adults$RealAge)),max(c(males_DNAmAge_adults$RealAge,females_DNAmAge_adults$RealAge)))
ylim=c(min(c(males_DNAmAge_adults$DNAmAge,females_DNAmAge_adults$DNAmAge)),max(c(males_DNAmAge_adults$DNAmAge,females_DNAmAge_adults$DNAmAge)))

plot(males_DNAmAge_adults$RealAge, males_DNAmAge_adults$DNAmAge, xlab='Chronological Age',ylab='DNAm Age', xlim=xlim, ylim=ylim,
  pch=21, col='black', bg='blue', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
abline(intercept,slope,lty='solid',col='blue',lwd=4)

fit=lm(females_DNAmAge_adults$DNAmAge~females_DNAmAge_adults$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(females_DNAmAge_adults$RealAge, females_DNAmAge_adults$DNAmAge, method='spearman', alternative='two.sided', exact=FALSE)
female_legend_label=paste0('Female: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

points(females_DNAmAge_adults$RealAge, females_DNAmAge_adults$DNAmAge,pch=21, col='black', bg='red', cex=2)
abline(intercept,slope,lty='solid',col='red',lwd=4)

abline(0,1,lty=2,col='black',lwd=4)

legend('topleft', c(male_legend_label, female_legend_label), cex=.8)
legend('bottomright', as.expression(bquote(bold(">21 y.o."))), cex=1.5, bty='n')

#plot the legend

# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("center", c('Males', 'Females','Identity'), lwd=4, col=c('blue','red','black'), lty=c(1,1,2), cex=1.3)

dev.off()


}

# Analysis of genes that lose DNAm in HGG

lost_sex_DMPs_interrogation = function(cases, controls, lost_DMPs_overlapping_genes_filename, 
  lost_DMPs_overlapping_genes_enrichR_filename, lost_DMPs_overlapping_imprinted_genes_filename, genes, ORA_table_filename) {

require(GenomicRanges)

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

mm=match(controls$probe, cases$probe)
loss=controls[which(is.na(mm)),]
loss=loss[,c(1,3,5,7:ncol(loss))]

rr=range(abs(loss$slope))
rr=paste(signif(rr[1],2),signif(rr[2],2),sep='-')

gr=GRanges(seqnames=loss$chr, 
  ranges=IRanges(loss$pos, loss$pos))

# Find genes that overlap or are <= 5kb from genes
oo=as.data.frame(findOverlaps(gr, genes, maxgap=5000, select='all', ignore.strand=TRUE))
hits=genes[oo$subjectHits,]
dist=distance(gr[oo$queryHits],hits,ignore.strand=TRUE)

pt1=loss[oo$queryHits,]

loss1=cbind(pt1[,-6], gene=as.vector(hits$Gene), Ensembl=as.vector(hits$Geneid), distance=dist)
rownames(loss1)=NULL

# How much of the lost DMPs overlap human TFs?
mm=match(loss1$Ensembl,TFs$Ensembl_ID)
temp=loss1[which(!is.na(mm)),]
mm=match(temp$Ensembl,TFs$Ensembl_ID)
num_TF_DMPs=length(unique(temp$probe))

# Regarding the queston above, how many unique TFs are there?
num_TFs=length(unique(mm))

loss1$TF=FALSE
mm=match(loss1$Ensembl,TFs$Ensembl_ID)
loss1$TF[!is.na(mm)]=TRUE

# How much of the lost DMPs overlap imprinted genes (genes that are imprinted in humans)?
mm=match(loss1$Ensembl,ig$Ensembl.ID)
temp=loss1[which(!is.na(mm)),]
mm=match(temp$Ensembl,ig$Ensembl.ID)
num_ig_DMPs=length(unique(temp$probe))

# Regarding the queston above, how many unique imprinted genes are there?
num_ig=length(unique(mm))

loss1$imprinted_gene=FALSE
mm=match(loss1$Ensembl,ig$Ensembl.ID)
loss1$imprinted_gene[!is.na(mm)]=TRUE

loss1$expressed_allele="NA"

index_ig=which(loss1$imprinted_gene==TRUE)

if (length(index_ig)>0){
  for (i in 1:length(index_ig)){
    m=match(loss1$Ensembl[index_ig[i]],ig$Ensembl.ID)
    loss1$expressed_allele[index_ig[i]]=as.vector(ig$ExpressedAllele[m])
  }
}

# Output the data

if(nrow(loss1) != 0){
  write.csv(loss1, file=lost_DMPs_overlapping_genes_filename, row.names=FALSE) # lost DMPs which overlap or are within 5kb of genes
  ug=data.frame(unique_genes=unique(as.vector(loss1$gene)))
  ug.Ensembl=unique(as.vector(loss1$Ensembl))
  write.csv(ug, file=lost_DMPs_overlapping_genes_enrichR_filename, row.names=FALSE)
}

# If any lost DMPs overlap an imprinted gene, make separate csv with the imprinted gene

if (num_ig_DMPs !=0 ){
  imp_gene_out=loss1[which(loss1$imprinted_gene=='TRUE'), c(4,5,2,3,8,9,10,13)]
  oo=order(imp_gene_out$fdr)
  imp_gene_out=imp_gene_out[oo,]
  write.csv(imp_gene_out, file=lost_DMPs_overlapping_imprinted_genes_filename, row.names=FALSE)
}

# run gene ontology on genes that overlap or are within 5kb of the 'lost sex-DMPs' ############################################
require("KEGG.db")
KEGG.PATHID2NAME=unlist(as.list(KEGGPATHID2NAME))
KEGG.PATHID2NAME=data.frame(kegg.path.ID=names(KEGG.PATHID2NAME), kegg.path.name=as.vector(KEGG.PATHID2NAME))

require(goseq)

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/gene_universe_450k.rda")
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/num_probes_per_gene.rda")

mm=match(universe.450k,ug.Ensembl)
indicator=rep(0, times=length(universe.450k))
indicator[which(!is.na(mm))]=1

aa=indicator
names(aa)=universe.450k

pwf=nullp(aa, "hg19", "ensGene", bias.data=num_probes_per_gene[,2])
GO.KEGG.wall=goseq(pwf,"hg19","ensGene",test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"))
GO.KEGG.wall$over_represented_FDR=p.adjust(GO.KEGG.wall$over_represented_pvalue, method="BH")
GO.KEGG.wall$ontology[which(is.na(GO.KEGG.wall$ontology))]='KEGG'

mm=match(as.vector(GO.KEGG.wall$category),as.vector(KEGG.PATHID2NAME$kegg.path.ID))

indexes=which(!is.na(mm))
mmm=mm[-which(is.na(mm))]

for (i in 1:length(indexes)) {
  GO.KEGG.wall$term[indexes[i]]=as.vector(KEGG.PATHID2NAME$kegg.path.name)[mmm[i]]
}

# from https://support.bioconductor.org/p/102273/
getGenes <- function(pwf, goterm, genome, ids){
    gene2cat <-  getgo(rownames(pwf), genome, ids, fetch.cats=c("GO:CC","GO:BP","GO:MF", "KEGG"))
    cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                      unlist(gene2cat, use.names = FALSE))
    out <- pwf[cat2gene[[goterm]],]
    out <- out[out$DEgenes > 0,]
    out
}

GO.KEGG.wall=GO.KEGG.wall[order(GO.KEGG.wall$over_represented_pvalue),]

nn =  length(which(GO.KEGG.wall$over_represented_pvalue<=0.05))

GO.KEGG.wall.out=GO.KEGG.wall[1:nn,]

GO.KEGG.wall.out$genes.symbols=NA
GO.KEGG.wall.out$genes.Ensembl=NA

for (i in 1:nrow(GO.KEGG.wall.out)){

  aa=getGenes(pwf, GO.KEGG.wall.out$category[i], "hg19", "ensGene")
  GO.KEGG.wall.out$genes.symbols[i]=paste(genes.df$Gene[match(rownames(aa), genes.df$Geneid)], collapse=", ")
  GO.KEGG.wall.out$genes.Ensembl[i]=paste(rownames(aa), collapse=", ")

  }

op=which(GO.KEGG.wall.out$ontology=='KEGG')

for(i in 1:length(op)){
  GO.KEGG.wall.out$category[op[i]]=paste0('KEGG ',GO.KEGG.wall.out$category[op[i]])
}

# GO.KEGG.wall.out=GO.KEGG.wall.out[,c(6,7,1,2,8,9,10)]

write.csv(GO.KEGG.wall.out, file=ORA_table_filename, row.names=FALSE)

#
dat=data.frame(
  num_DMPs_lost=nrow(loss),
  proportion_DMPs_lost=nrow(loss)/nrow(controls),
  range_of_abs_value_of_DNAm_differenes_between_sexes=rr,
  how_many_lost_DMPs_overlap_genes_or_promoters_or_enhancers= nrow(loss[which(loss$distance_to_nearest_gene<=5000 | loss$Enhancer == 'TRUE' | loss$Promoter == 'TRUE'),]),
  how_many_lost_DMPs_overlap_TFs=num_TF_DMPs,
  how_many_TFs=num_TFs,
  how_many_lost_DMPs_overlap_imprinted_genes=num_ig_DMPs,
  how_many_imprinted_genes=num_ig,
  number_sig_GO_or_KEGG_terms_FDR_0.1=length(which(GO.KEGG.wall$over_represented_FDR<=0.1)),
  number_sig_GO_or_KEGG_terms_unadjP_0.05=length(which(GO.KEGG.wall$over_represented_pvalue<=0.05)) 
)

dat=t(dat)

return(dat)

}

# modified minfi functions

getSnpInfo_mine = function (object, snpAnno = NULL)
{
    av <- minfi:::.availableAnnotation(object)
    if (is.null(snpAnno)) {
        snpAnno <- grep(pattern = "^SNPs\\.", x = getAnnotationObject(object)@defaults,
            value = TRUE)
    } else {
        snpAnno <- sub("^SNPs\\.", "", snpAnno)
        if (!snpAnno %in% av$classChoices$SNPs) {
            stop(sprintf("snpAnno '%s' is not part of the annotation",
                snpAnno))
        } else {
            snpAnno <- sprintf("SNPs.%s", snpAnno)
        }
    }
    snps <- getAnnotation(object, what = snpAnno)
    snps
}

dropLociWithSnps_mine= function (object, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)
{
    minfi:::.isGenomicOrStop(object)
    maf_cols <- paste0(snps, "_maf")
    snpDF <- getSnpInfo_mine(object, snpAnno = snpAnno)
    choices <- c("Probe_maf", "CpG_maf", "SBE_maf")
    if (!all(choices %in% colnames(snpDF))) {
        stop("The specificed 'snpAnno' is not supported by this function")
    }
    if (sum(!(maf_cols %in% choices)) > 0) {
        stop("snps vector argument must be a combination of  \"Probe\", ",
            "\"CpG\" and \"SBE\"")
    }
    if (!is.numeric(maf) || maf < 0 || maf > 1) {
        stop("maf argument must be a numeric value between 0 and 1")
    }
    wh <- Reduce(union, lapply(maf_cols, function(xx) {
        which(snpDF[, xx] >= maf)
    }))
    wh <- sort(wh)
    if (length(wh) == 0)
        return(object)
    object[-wh, ]
}


############################
## DRAFTS & ARCHIVE
############################

# ## plot blocks

# ## cset: result of minfi::cpgCollapse()$object
# ## blocks450: results of blockFinder
# ## coi = covariate of interest
# ## N: the number of blocks to plot, default=10
# ## blockname = name of block track, default='coi'
# ## filename = where to save plots
# ## scale, in kb. default = 100

# blockPlot = function(cset,
#   blocks450,
#   coi,
#   N=100,
#   blockname = "coi",
#   filename,
#   scale=100,
#   showMethPanel = TRUE,
#   showGenePanel=TRUE,
#   showDiffPanel=TRUE, 
#   showCancerPanel = FALSE,
#   bty= "o") {

#   panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
  
#   #require(GenomicRanges)

#   blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
#   colIds=match(c("chr","start","end"),names(blocks450$table))
#   mcols(blocksTable)=blocks450$table[-colIds]

#   plotRegion = blocksTable[1:N]

#   ## annotation based on ensembl
#   cat("Loading Annotation.\n")
#   #load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/genomicState.rda")
#   load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda') # txdb.gencode.v35.hg19, genomicState.gencode.v35.hg19, genes.gencode.v35.hg19, genes.df.gencode.v35.hg19
#   genomicState = genomicState.gencode.v35.hg19
#   gs = genomicState$fullGenome
#   oo = findOverlaps(blocksTable, gs)
#   anno = split(gs[subjectHits(oo)], queryHits(oo))

#   # cancer blocks
#   if(showCancerPanel) {
#     # read in cancer DNAm bocks (Hansen et al, 2011, PMID: 21706001)
#     #load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
#     blocks=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/cancer_blocks_hansen.csv')
#     blocks=GRanges(blocks)
#     genome(blocks)="hg19"
#     cancerBlocks = blocks
#   }

#   cat("Ploting.")
#   pdf(filename,height=5,width=10)
#   par(bty=bty)

#   for(i in seq(along=plotRegion)) { #//
#     cat(".")
#     r = plotRegion[i]
#     tmp=subsetByOverlaps(cset,r)
#     tmp450=sort(subsetByOverlaps(blocksTable,r))
#     if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
                      
#     beta=getBeta(tmp)
#     x=start(tmp)

#     ii=cset %over% r
#     d=blocks450$coef[ii]
#     sd=blocks450$fitted[ii]
#     ## which rows
#     Index=split(seq_along(coi),coi)
#     mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
#     smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

#     ## paneling, from hector
#     mypar(1,1, brewer.name = "Set1")
#     par(mar=par()$mar+c(0,3,0,0))
#     omar=par()$mar
#     cmar=omar
#     cmar[1]=.5
#     par(mar=cmar)
    
#     layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
#     if(showMethPanel) {
#       matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
#       xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
#       matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
#       legend("bottomright",col=seq(along=levels(factor(coi))),
#         lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
        
#       segments(min(x), .2, min(x)+scale*1000, .2)
#       text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
#       axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
#       mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
#       mtext("Methylation",side=2, line = 2.5,cex=1.5)

#       legend("topright", paste0("fwer = ",
#         signif(blocks450$tab$fwer[i],3)),bty=bty)
#     }   
#     # annotation
#     if(showGenePanel) {
#       cmar=omar
#       if(!is.na(showDiffPanel)) {
#         cmar[3]=0.5
#       } else {
#         cmar[c(1,3)]=c(0.5,0.5)
#       }
#       par(mar=cmar)

#       plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
#          xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
#       a = as.data.frame(anno[[i]])
#       Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
#       Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
#       Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
#       axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
#       abline(h=0,lty=3)
#       for(k in 1:nrow(a)) {
#         polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
#           Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
#       }

#       ## by gene
#       g = split(a, sapply(a$symbol,"[", 1))
#       # g = split(a, a$Symbol)
#       s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
#       g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
      
#       if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
          
#       mtext("Genes",side=2, line = 2.5,cex=1.5)
#       if(!showDiffPanel) {
#         xtick=pretty(x)
#         axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
#       }
#     }

#     ## mean diff
#     if(showDiffPanel) {
#       cmar=omar
#       cmar[c(1,3)]=c(.5,.5)
#       par(mar=cmar)

#       zz=granges(tmp)

#       matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
#       axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
#       xtick=pretty(x)
#       axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
#       mtext("Diff",side=2, line = 2.5,cex=1.5)

#       ii=which(zz$type=="OpenSea")
#       blockgroup=zz$blockgroup[ii]

#       blockIndexes=split(seq(along=blockgroup),blockgroup)
#       for (ind in blockIndexes) {
#         ind=ii[ind]
#         lines(x[ind], sd[ind], lwd=2.5,col="black")
#       }


#       points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
#       axis(side=2,at=c(-2,0,2))
#       abline(h=0,lty=2,col="black")

#       cmar=omar
#       cmar[3]=.5
#       par(mar=cmar)
#       matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
#         ylab="",ylim=c(0,2),bty="n")
#     }
    
#     #  browser()
#     if(showCancerPanel) {
#       col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
#       rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
#       if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$Direction.of.Methylation.Change=="hypo","blue","red"))
#       axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
#       legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
#     }
#   } #//
  
#   dev.off()

# }

# ## cset: result of cpgCollapse()$cobj
# ## block450: results of blockFinder
# ## coi = covariate of interest
# ## N: the number of blocks to plot, default=10
# ## blockname = name of block track, default='coi'
# ## filename = where to save plots
# ## scale, in kb. default = 100
# blockPlot = function(cset, blocks450, coi, N=10,
#   blockname = "coi", filename=paste0(blockname,"_blocks.pdf"),scale=100,
#   showMethPanel = TRUE, showGenePanel=TRUE, showDiffPanel=TRUE, 
#   showCancerPanel = TRUE, bty= "o") {
#   panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
  
#   require(GenomicRanges)

#   blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
#   colIds=match(c("chr","start","end"),names(blocks450$table))
#   mcols(blocksTable)=blocks450$table[-colIds]

#   plotRegion = blocksTable[1:N]

#   ## annotation based on ensembl
#   cat("Loading Annotation.\n")
#   load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
#   gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
#   oo = findOverlaps(blocksTable, gs)
#   anno = split(gs[subjectHits(oo)], queryHits(oo))

#   # cancer blocks
#   if(showCancerPanel) {
#     load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
#     genome(blocks)="hg19"
#     cancerBlocks = blocks
#   }

#   cat("Ploting.")
#   pdf(filename,height=5,width=10)
#   par(bty=bty)

#   for(i in seq(along=plotRegion)) {
#     cat(".")
#     r = plotRegion[i]
#     tmp=subsetByOverlaps(cset,r)
#     tmp450=sort(subsetByOverlaps(blocksTable,r))
#     if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
                      
#     beta=getBeta(tmp)
#     x=start(tmp)

#     ii=cset %over% r
#     d=blocks450$coef[ii]
#     sd=blocks450$fitted[ii]
#     ## which rows
#     Index=split(seq_along(coi),coi)
#     mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
#     smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

#     ## paneling, from hector
#     mypar(1,1, brewer.name = "Set1")
#     par(mar=par()$mar+c(0,3,0,0))
#     omar=par()$mar
#     cmar=omar
#     cmar[1]=.5
#     par(mar=cmar)
    
#     layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
#     if(showMethPanel) {
#       matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
#       xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
#       matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
#       legend("bottomright",col=seq(along=levels(factor(coi))),
#         lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
        
#       segments(min(x), .2, min(x)+scale*1000, .2)
#       text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
#       axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
#       mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
#       mtext("Methylation",side=2, line = 2.5,cex=1.5)

#       legend("topright", paste0("fwer = ",
#         signif(blocks450$tab$fwer[i],3)),bty=bty)
#     }   
#     # annotation
#     if(showGenePanel) {
#       cmar=omar
#       if(!is.na(showDiffPanel)) {
#         cmar[3]=0.5
#       } else {
#         cmar[c(1,3)]=c(0.5,0.5)
#       }
#       par(mar=cmar)

#       plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
#          xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
#       a = as.data.frame(anno[[i]])
#       Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
#       Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
#       Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
#       axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
#       abline(h=0,lty=3)
#       for(k in 1:nrow(a)) {
#         polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
#           Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
#       }

#       ## by gene
#       g = split(a, sapply(a$symbol,"[", 1))
#       # g = split(a, a$Symbol)
#       s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
#       g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
      
#       if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
          
#       mtext("Genes",side=2, line = 2.5,cex=1.5)
#       if(!showDiffPanel) {
#         xtick=pretty(x)
#         axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
#       }
#     }

#     ## mean diff
#     if(showDiffPanel) {
#       cmar=omar
#       cmar[c(1,3)]=c(.5,.5)
#       par(mar=cmar)

#       zz=granges(tmp)

#       matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
#       axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
#       xtick=pretty(x)
#       axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
#       mtext("Diff",side=2, line = 2.5,cex=1.5)

#       ii=which(zz$type=="OpenSea")
#       blockgroup=zz$blockgroup[ii]

#       blockIndexes=split(seq(along=blockgroup),blockgroup)
#       for (ind in blockIndexes) {
#         ind=ii[ind]
#         lines(x[ind], sd[ind], lwd=2.5,col="black")
#       }


#       points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
#       axis(side=2,at=c(-2,0,2))
#       abline(h=0,lty=2,col="black")

#       cmar=omar
#       cmar[3]=.5
#       par(mar=cmar)
#       matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
#         ylab="",ylim=c(0,2),bty="n")
#     }
    
#     #  browser()
#     if(showCancerPanel) {
#       col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
#       rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
#       if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$direction=="hypo","blue","red"))
#       axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
#       legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
#     }
#   }
#   dev.off()
# }
