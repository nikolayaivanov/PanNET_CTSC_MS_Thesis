#
rm(list=ls())
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

## macrophage data
macrophage_frequencies=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/Murine_FACS/macrophage_analysis.csv', header=TRUE)
macrophage_frequencies$Live.CD45.F480.CD206hi=macrophage_frequencies$Live.CD45.F480.Q2+macrophage_frequencies$Live.CD45.F480.Q3
macrophage_frequencies$Live.CD45.F480.CD206lo=macrophage_frequencies$Live.CD45.F480.Q1+macrophage_frequencies$Live.CD45.F480.Q4

# normalize by frequency of Live/CD45+ cells
#normalization_factor=macrophage_frequencies$Live.CD45*macrophage_frequencies$Live.CD45.F480
normalization_factor=macrophage_frequencies$Live.CD45.F480
for (i in 4:ncol(macrophage_frequencies)){
	macrophage_frequencies[,i]=macrophage_frequencies[,i]/normalization_factor
}

## T-cell data
T_cell_frequencies=read.csv('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/Murine_FACS/T_cell_analysis.csv', header=TRUE)

# normalize by frequency of Live/CD45+ cells
#normalization_factor=T_cell_frequencies$Live.CD45*T_cell_frequencies$Live.CD45.CD3
#normalization_factor=T_cell_frequencies$Live.CD45.CD3
T_cell_frequencies$Live.CD45.CD3.CD4 = T_cell_frequencies$Live.CD45.CD3.CD4/T_cell_frequencies$Live.CD45.CD3
T_cell_frequencies$Live.CD45.CD3.CD4.CTLA4 = T_cell_frequencies$Live.CD45.CD3.CD4.CTLA4/T_cell_frequencies$Live.CD45.CD3.CD4
T_cell_frequencies$Live.CD45.CD3.CD8 = T_cell_frequencies$Live.CD45.CD3.CD8/T_cell_frequencies$Live.CD45.CD3
T_cell_frequencies$Live.CD45.CD3.CD8.CTLA4 = T_cell_frequencies$Live.CD45.CD3.CD8.CTLA4/T_cell_frequencies$Live.CD45.CD3.CD8

# for (i in 4:ncol(T_cell_frequencies)){
# 	T_cell_frequencies[,i]=T_cell_frequencies[,i]/normalization_factor
# }

#pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/Murine_FACS/pdf/FACS_plots_full_page_normalized_by_Live.CD45.F480_OR_Live.CD45.CD3_AND_Live.CD45.pdf')
pdf('/athena/masonlab/scratch/users/nai2008/PNET_thesisProject/Murine_FACS/pdf/FACS_plots_full_page_normalized.pdf')

###################################################################################
## T-cells
###################################################################################

# pancreas indicies
index_riptag_PNET=grep('riptag_PNET',T_cell_frequencies$Sample)
	riptag_PNET_mouseID = ss(T_cell_frequencies$Sample[index_riptag_PNET],'_',1)

index_riptag_pancreas=grep('riptag_pancreas',T_cell_frequencies$Sample)
	riptag_pancreas_mouseID = ss(T_cell_frequencies$Sample[index_riptag_pancreas],'_',1)

mm=match(riptag_PNET_mouseID,riptag_pancreas_mouseID)
if(length(which(is.na(mm))) != 0 ){ 
	index_riptag_PNET=index_riptag_PNET[-which(is.na(mm))]
	riptag_PNET_mouseID = ss(T_cell_frequencies$Sample[index_riptag_PNET],'_',1)
}

mm=match(riptag_PNET_mouseID,riptag_pancreas_mouseID)
index_riptag_pancreas=index_riptag_pancreas[mm]

riptag_PNET_mouseID = ss(T_cell_frequencies$Sample[index_riptag_PNET],'_',1)
riptag_pancreas_mouseID = ss(T_cell_frequencies$Sample[index_riptag_pancreas],'_',1)
all(riptag_PNET_mouseID == riptag_pancreas_mouseID) #TRUE

index_control_pancreas=grep('control_pancreas',T_cell_frequencies$Sample)

# blood indicies
index_riptag_blood=grep('riptag_blood',T_cell_frequencies$Sample)
index_control_blood=grep('control_blood',T_cell_frequencies$Sample)

# spleen indicies
index_riptag_spleen=grep('riptag_spleen',T_cell_frequencies$Sample)
index_control_spleen=grep('control_spleen',T_cell_frequencies$Sample)

par(mar=c(5,5.3,4.1,2.1))

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("T-cell Flow Cytometry Data\n",'(Mouse models)'), cex = 1.6, col = "black")

####################
## Live/CD45+ Cells
####################

### pancreas
riptag_pancreas=T_cell_frequencies$Live.CD45[index_riptag_pancreas]
riptag_PNET=T_cell_frequencies$Live.CD45[index_riptag_PNET]
control_pancreas=T_cell_frequencies$Live.CD45[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Cell frequency (%)',
	names=c('PNET','RipTag pancreas','Control pancreas'), main=c('Live/CD45+ Cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=T_cell_frequencies$Live.CD45[index_riptag_blood]
control_blood=T_cell_frequencies$Live.CD45[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Cell frequency (%)',names=c('RipTag blood','Control blood'),main=c('Live/CD45+ Cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=T_cell_frequencies$Live.CD45[index_riptag_spleen]
control_spleen=T_cell_frequencies$Live.CD45[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Cell frequency (%)',names=c('RipTag spleen','Control spleen'),main=c('Live/CD45+ Cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## Live/CD45+/CD3+ Cells
####################

### pancreas
riptag_pancreas=T_cell_frequencies$Live.CD45.CD3[index_riptag_pancreas]
riptag_PNET=T_cell_frequencies$Live.CD45.CD3[index_riptag_PNET]
control_pancreas=T_cell_frequencies$Live.CD45.CD3[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Cell frequency (%)',
	names=c('PNET','RipTag pancreas','Control pancreas'), main=c('Live/CD45+/CD3+ Cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=T_cell_frequencies$Live.CD45.CD3[index_riptag_blood]
control_blood=T_cell_frequencies$Live.CD45.CD3[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Cell frequency (%)',names=c('RipTag blood','Control blood'),main=c('Live/CD45+/CD3+ Cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=T_cell_frequencies$Live.CD45.CD3[index_riptag_spleen]
control_spleen=T_cell_frequencies$Live.CD45.CD3[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Cell frequency (%)',names=c('RipTag spleen','Control spleen'),main=c('Live/CD45+/CD3+ Cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## CD8+ T-cells
####################

### pancreas
riptag_pancreas=T_cell_frequencies$Live.CD45.CD3.CD8[index_riptag_pancreas]
riptag_PNET=T_cell_frequencies$Live.CD45.CD3.CD8[index_riptag_PNET]
control_pancreas=T_cell_frequencies$Live.CD45.CD3.CD8[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',
	names=c('PNET','RipTag pancreas','Control pancreas'), main=c('CD8+ T-cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=T_cell_frequencies$Live.CD45.CD3.CD8[index_riptag_blood]
control_blood=T_cell_frequencies$Live.CD45.CD3.CD8[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('CD8+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=T_cell_frequencies$Live.CD45.CD3.CD8[index_riptag_spleen]
control_spleen=T_cell_frequencies$Live.CD45.CD3.CD8[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('CD8+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## CD4+ T-cells
####################

### pancreas
riptag_pancreas=T_cell_frequencies$Live.CD45.CD3.CD4[index_riptag_pancreas]
riptag_PNET=T_cell_frequencies$Live.CD45.CD3.CD4[index_riptag_PNET]
control_pancreas=T_cell_frequencies$Live.CD45.CD3.CD4[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('CD4+ T-cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=T_cell_frequencies$Live.CD45.CD3.CD4[index_riptag_blood]
control_blood=T_cell_frequencies$Live.CD45.CD3.CD4[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('CD4+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=T_cell_frequencies$Live.CD45.CD3.CD4[index_riptag_spleen]
control_spleen=T_cell_frequencies$Live.CD45.CD3.CD4[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('CD4+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## CTLA+ CD8+ T-cells
####################

### pancreas
riptag_pancreas=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_riptag_pancreas]
riptag_PNET=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_riptag_PNET]
control_pancreas=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('CD8+/CTLA4+ T-cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

boxplot(riptag_PNET,riptag_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','Non-neoplastic Pancreas'), main=c('CD8+/CTLA4+ T-cells'),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
legend(x='topleft',legend=paste0('p = ', p_riptag_PNET_vs_riptag_pancreas), bty='n')
panc.all.dat=panc.all.dat[1:length(c(riptag_PNET,riptag_pancreas))]
panc.all.dat.labs=panc.all.dat.labs[1:length(c(riptag_PNET,riptag_pancreas))]
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_riptag_blood]
control_blood=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('CD8+/CTLA4+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_riptag_spleen]
control_spleen=T_cell_frequencies$Live.CD45.CD3.CD8.CTLA[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('CD8+/CTLA4+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## CTLA+ CD4+ T-cells
####################

### pancreas
riptag_pancreas=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_riptag_pancreas]
riptag_PNET=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_riptag_PNET]
control_pancreas=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('CD4+/CTLA4+ T-cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

boxplot(riptag_PNET,riptag_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','Non-neoplastic pancreas'), main=c('CD4+/CTLA4+ T-cells'),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
legend(x='topleft',legend=paste0('p = ', p_riptag_PNET_vs_riptag_pancreas), bty='n')
panc.all.dat=panc.all.dat[1:length(c(riptag_PNET,riptag_pancreas))]
panc.all.dat.labs=panc.all.dat.labs[1:length(c(riptag_PNET,riptag_pancreas))]
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

boxplot(control_pancreas, riptag_PNET, col='grey93', ylab='Normalized cell frequency',names=c('Control pancreas', 'PNET'), main=c('CD4+/CTLA4+ T-cells'),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
legend(x='topleft',legend=paste0('p = ', p_riptag_PNET_vs_control_pancreas), bty='n')
panc.all.dat=c(control_pancreas, riptag_PNET)
panc.all.dat.labs=c(rep('control_pancreas',times=length(control_pancreas)),rep('riptag_PNET', times=length(riptag_PNET)))
panc.all.dat.labs=factor(panc.all.dat.labs,levels=c('control_pancreas','riptag_PNET'))
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_riptag_blood]
control_blood=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('CD4+/CTLA4+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_riptag_spleen]
control_spleen=T_cell_frequencies$Live.CD45.CD3.CD4.CTLA[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('CD4+/CTLA4+ T-cells',p),cex.main=1,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)


###################################################################################
## Macrophages
###################################################################################

# pancreas indicies
index_riptag_PNET=grep('riptag_PNET',macrophage_frequencies$Sample)
	riptag_PNET_mouseID = ss(macrophage_frequencies$Sample[index_riptag_PNET],'_',1)

index_riptag_pancreas=grep('riptag_pancreas',macrophage_frequencies$Sample)
	riptag_pancreas_mouseID = ss(macrophage_frequencies$Sample[index_riptag_pancreas],'_',1)

mm=match(riptag_PNET_mouseID,riptag_pancreas_mouseID)
if(length(which(is.na(mm))) != 0 ){ 
	index_riptag_PNET=index_riptag_PNET[-which(is.na(mm))]
	riptag_PNET_mouseID = ss(macrophage_frequencies$Sample[index_riptag_PNET],'_',1)
}

mm=match(riptag_PNET_mouseID,riptag_pancreas_mouseID)
index_riptag_pancreas=index_riptag_pancreas[mm]

riptag_PNET_mouseID = ss(macrophage_frequencies$Sample[index_riptag_PNET],'_',1)
riptag_pancreas_mouseID = ss(macrophage_frequencies$Sample[index_riptag_pancreas],'_',1)
all(riptag_PNET_mouseID == riptag_pancreas_mouseID) #TRUE

index_control_pancreas=grep('control_pancreas',macrophage_frequencies$Sample)

# blood indicies
index_riptag_blood=grep('riptag_blood',macrophage_frequencies$Sample)
index_control_blood=grep('control_blood',macrophage_frequencies$Sample)

# spleen indicies
index_riptag_spleen=grep('riptag_spleen',macrophage_frequencies$Sample)
index_control_spleen=grep('control_spleen',macrophage_frequencies$Sample)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("Macrophage Flow Cytometry Data\n",'(Mouse models)'), cex = 1.6, col = "black")

####################
## Live/CD45+ Cells
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Cell frequency (%)',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('Live/CD45+ Cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Cell frequency (%)',names=c('RipTag blood','Control blood'),main=c('Live/CD45+ Cells',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Cell frequency (%)',names=c('RipTag spleen','Control spleen'),main=c('Live/CD45+ Cells',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## Live/CD45+/F480+ Cells
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45.F480[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45.F480[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45.F480[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Cell frequency (%)',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('Live/CD45+/F480+ Cells',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45.F480[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45.F480[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Cell frequency (%)',names=c('RipTag blood','Control blood'),main=c('Live/CD45+/F480+ Cells',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45.F480[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45.F480[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Cell frequency (%)',names=c('RipTag spleen','Control spleen'),main=c('Live/CD45+/F480+ Cells',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## M1 Macrophages (Q1) [CD11C+/CD206-]
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45.F480.Q1[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45.F480.Q1[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45.F480.Q1[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('M1 Macrophages (Q1) [CD11C+/CD206-]',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45.F480.Q1[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45.F480.Q1[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('M1 Macrophages (Q1) [CD11C+/CD206-]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45.F480.Q1[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45.F480.Q1[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('M1 Macrophages (Q1) [CD11C+/CD206-]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## Macrophages (Q2) [CD11C+/CD206+]
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45.F480.Q2[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45.F480.Q2[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45.F480.Q2[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('Macrophages (Q2) [CD11C+/CD206+]',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45.F480.Q2[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45.F480.Q2[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('Macrophages (Q2) [CD11C+/CD206+]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45.F480.Q2[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45.F480.Q2[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('Macrophages (Q2) [CD11C+/CD206+]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## M2 Macrophages (Q3) [CD11C-/CD206+]
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45.F480.Q3[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45.F480.Q3[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45.F480.Q3[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('M2 Macrophages (Q3) [CD11C-/CD206+]',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45.F480.Q3[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45.F480.Q3[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('M2 Macrophages (Q3) [CD11C-/CD206+]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45.F480.Q3[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45.F480.Q3[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('M2 Macrophages (Q3) [CD11C-/CD206+]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## Macrophages (Q4) [CD11C-/CD206-]
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45.F480.Q4[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45.F480.Q4[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45.F480.Q4[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('Macrophages (Q4) [CD11C-/CD206-]',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45.F480.Q4[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45.F480.Q4[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('Macrophages (Q4) [CD11C-/CD206-]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45.F480.Q4[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45.F480.Q4[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('Macrophages (Q4) [CD11C-/CD206-]',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

####################
## CD206 hi macrophages
####################

#pancreas
riptag_pancreas=macrophage_frequencies$Live.CD45.F480.CD206hi[index_riptag_pancreas]
riptag_PNET=macrophage_frequencies$Live.CD45.F480.CD206hi[index_riptag_PNET]
control_pancreas=macrophage_frequencies$Live.CD45.F480.CD206hi[index_control_pancreas]

panc.all.dat=c(riptag_PNET,riptag_pancreas,control_pancreas)
panc.all.dat.labs=c(rep('riptag_PNET',times=length(riptag_PNET)),rep('riptag_pancreas',times=length(riptag_pancreas)),rep('control_pancreas',times=length(control_pancreas)))
panc.all.dat.labs=factor(panc.all.dat.labs, levels=c('riptag_PNET','riptag_pancreas','control_pancreas'))

# calculate RipTag pancreas vs control pancreas p-value (Wilcoxon rank-sum (Mann-Whitney))
#p_riptag_pancreas_vs_control_pancreas=signif(wilcox.test(riptag_pancreas, control_pancreas)$p.value,2)
# t-test:
#p_riptag_pancreas_vs_control_pancreas=signif(t.test(riptag_pancreas,control_pancreas)$p.value,2)

# calculate PNET vs control pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_control_pancreas=signif(wilcox.test(riptag_PNET,control_pancreas, paired=TRUE)$p.value,2)
# t-test:
p_riptag_PNET_vs_control_pancreas=signif(t.test(riptag_PNET,control_pancreas)$p.value,2)

# caluclate PNET vs RipTag pancreas p-value (Wilcoxon signed-rank)
#p_riptag_PNET_vs_riptag_pancreas=signif(wilcox.test(riptag_PNET, riptag_pancreas, paired=TRUE)$p.value,2)
# paired t-test:
p_riptag_PNET_vs_riptag_pancreas=signif(t.test(riptag_PNET,riptag_pancreas, paired=TRUE)$p.value,2)

FDR=signif(p.adjust(c(p_riptag_PNET_vs_control_pancreas, p_riptag_PNET_vs_riptag_pancreas),method='BH'),2)

zag2=paste0('PNET vs Control pancreas: p = ',p_riptag_PNET_vs_control_pancreas,'; FDR = ',FDR[1])
zag3=paste0('PNET vs matched non-neoplastic pancreas: p = ',p_riptag_PNET_vs_riptag_pancreas,'; FDR = ',FDR[2])
boxplot(riptag_PNET,riptag_pancreas,control_pancreas, col='grey93', ylab='Normalized cell frequency',names=c('PNET','RipTag pancreas','Control pancreas'), main=c('CD206 high macrophages',zag2,zag3),cex.main=.6,las=1,outline=FALSE,cex.lab=1.5)
points(panc.all.dat ~ jitter(as.numeric(panc.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#blood
riptag_blood=macrophage_frequencies$Live.CD45.F480.CD206hi[index_riptag_blood]
control_blood=macrophage_frequencies$Live.CD45.F480.CD206hi[index_control_blood]

blood.all.dat=c(riptag_blood,control_blood)
blood.all.dat.labs=c(rep('riptag_blood',times=length(riptag_blood)),rep('control_blood',times=length(control_blood)))
blood.all.dat.labs=factor(blood.all.dat.labs, levels=c('riptag_blood','control_blood'))

#p=paste('p =',signif(wilcox.test(riptag_blood,control_blood)$p.value,2))
p=paste('p =',signif(t.test(riptag_blood,control_blood)$p.value,2))
boxplot(riptag_blood,control_blood, col='grey93', ylab='Normalized cell frequency',names=c('RipTag blood','Control blood'),main=c('CD206 high macrophages',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(blood.all.dat ~ jitter(as.numeric(blood.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

#spleen
riptag_spleen=macrophage_frequencies$Live.CD45.F480.CD206hi[index_riptag_spleen]
control_spleen=macrophage_frequencies$Live.CD45.F480.CD206hi[index_control_spleen]

spleen.all.dat=c(riptag_spleen,control_spleen)
spleen.all.dat.labs=c(rep('riptag_spleen',times=length(riptag_spleen)),rep('control_spleen',times=length(control_spleen)))
spleen.all.dat.labs=factor(spleen.all.dat.labs, levels=c('riptag_spleen','control_spleen'))

#p=paste('p =',signif(wilcox.test(riptag_spleen,control_spleen)$p.value,2))
p=paste('p =',signif(t.test(riptag_spleen,control_spleen)$p.value,2))
boxplot(riptag_spleen,control_spleen, col='grey93', ylab='Normalized cell frequency',names=c('RipTag spleen','Control spleen'),main=c('CD206 high macrophages',p),cex.main=.9,las=1,outline=FALSE,cex.lab=1.5)
points(spleen.all.dat ~ jitter(as.numeric(spleen.all.dat.labs),amount=0.1), pch=21, col='black', bg='darkgrey',cex=1.2)

dev.off()

# NAIvanov

