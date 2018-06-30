pol2_0_TSS = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TSS_ALL_-3000_+0\\RN6_TSS_ALL_-3000_+0_RASMC_POL2_UNSTIM_NEW.gff',header=TRUE)
pol2_0_TXN = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TXN_ALL_-0_+0\\RN6_TXN_ALL_-0_+0_RASMC_POL2_UNSTIM_NEW.gff',header=TRUE)
pol2_0_TTR = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TTR_ALL_-0_+3000\\RN6_TTR_ALL_-0_+3000_RASMC_POL2_UNSTIM_NEW.gff',header=TRUE)

pol2_2_TSS = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TSS_ALL_-3000_+0\\RN6_TSS_ALL_-3000_+0_RASMC_POL2_PDGF_2H_NEW.gff',header=TRUE)
pol2_2_TXN = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TXN_ALL_-0_+0\\RN6_TXN_ALL_-0_+0_RASMC_POL2_PDGF_2H_NEW.gff',header=TRUE)
pol2_2_TTR = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TTR_ALL_-0_+3000\\RN6_TTR_ALL_-0_+3000_RASMC_POL2_PDGF_2H_NEW.gff',header=TRUE)

pol2_24_TSS = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TSS_ALL_-3000_+0\\RN6_TSS_ALL_-3000_+0_RASMC_POL2_PDGF_24H_NEW.gff',header=TRUE)
pol2_24_TXN = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TXN_ALL_-0_+0\\RN6_TXN_ALL_-0_+0_RASMC_POL2_PDGF_24H_NEW.gff',header=TRUE)
pol2_24_TTR = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TTR_ALL_-0_+3000\\RN6_TTR_ALL_-0_+3000_RASMC_POL2_PDGF_24H_NEW.gff',header=TRUE)

pol2_2J_TSS = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TSS_ALL_-3000_+0\\RN6_TSS_ALL_-3000_+0_RASMC_POL2_PDGF_2H_JQ1_NEW.gff',header=TRUE)
pol2_2J_TXN = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TXN_ALL_-0_+0\\RN6_TXN_ALL_-0_+0_RASMC_POL2_PDGF_2H_JQ1_NEW.gff',header=TRUE)
pol2_2J_TTR = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TTR_ALL_-0_+3000\\RN6_TTR_ALL_-0_+3000_RASMC_POL2_PDGF_2H_JQ1_NEW.gff',header=TRUE)

pol2_24J_TSS = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TSS_ALL_-3000_+0\\RN6_TSS_ALL_-3000_+0_RASMC_POL2_PDGF_24H_JQ1_NEW.gff',header=TRUE)
pol2_24J_TXN = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TXN_ALL_-0_+0\\RN6_TXN_ALL_-0_+0_RASMC_POL2_PDGF_24H_JQ1_NEW.gff',header=TRUE)
pol2_24J_TTR = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\mappedFolder\\RN6_TTR_ALL_-0_+3000\\RN6_TTR_ALL_-0_+3000_RASMC_POL2_PDGF_24H_JQ1_NEW.gff',header=TRUE)

active_genes = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\activeListTable.txt', header = FALSE)

orderTable = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\170206_waterfall_genes_log2_24_v_0_ordered_BRD4_H3k_RNA_0_2_24_and_log2.txt',header = TRUE)
up_24=which(orderTable[,12]>=1)
down_24 = which(orderTable[,12]<=-1)

orderTable2 = read.delim('C:\\Users\\rhirsch\\Documents\\rasmc_docs\\170206_waterfall_genes_log2_2_v_0_ordered_BRD4_H3k_RNA_0_2_24_and_log2.txt',header = TRUE)
up_2=which(orderTable2[,10]>=1)
down_2 = which(orderTable2[,10]<=-1)


pol2_0 = cbind(pol2_0_TSS[,3:62],pol2_0_TXN[,3:202],pol2_0_TTR[,3:62])
rownames(pol2_0) = pol2_0_TTR[,1]

pol2_2 = cbind(pol2_2_TSS[,3:62], pol2_2_TXN[,3:202],pol2_2_TTR[,3:62])
rownames(pol2_2) = rownames(pol2_0)

pol2_24 = cbind(pol2_24_TSS[,3:62],pol2_24_TXN[,3:202],pol2_24_TTR[,3:62])
rownames(pol2_24) = pol2_0_TTR[,1]

pol2_2J = cbind(pol2_2J_TSS[,3:62], pol2_2J_TXN[,3:202],pol2_2J_TTR[,3:62])
rownames(pol2_2J) = rownames(pol2_0)

pol2_24J = cbind(pol2_24J_TSS[,3:62], pol2_24J_TXN[,3:202],pol2_24J_TTR[,3:62])
rownames(pol2_24J) = rownames(pol2_0)


###########################################
############Active Genes Metas#############
###########################################

active0 = c()
for(i in 1:length(active_genes[,1])){
  row = which(as.character(rownames(pol2_0))==as.character(active_genes[i,1]))
  active0=c(active0,row)
}

pol2_0_meta = apply(pol2_0[active0,],2,mean,na.rm=TRUE)
pol2_2_meta = apply(pol2_2[active0,],2,mean,na.rm=TRUE)
pol2_24_meta = apply(pol2_24[active0,],2,mean,na.rm=TRUE)
pol2_2J_meta = apply(pol2_2J[active0,],2,mean,na.rm=TRUE)
pol2_24J_meta = apply(pol2_24J[active0,],2,mean,na.rm=TRUE)


pdf(file='C:\\Users\\rhirsch\\Documents\\rasmc_figures_1-3\\170222_rasmc_active_pol2_0v2_meta.pdf',width = 10,height = 8)
plot(1:320, pol2_0_meta,type='l',col='red',ylim =c(0,2.5),xaxt='n',xlab='',ylab='rpm/bp',main='active genes 0H (red) vs. 2H (black)')
lines(1:320, pol2_2_meta,type='l',col='black')
axis(1,c(0,60,260,320),c('-3kb','TSS','END','+3kb'))
dev.off()

pdf(file='C:\\Users\\rhirsch\\Documents\\rasmc_figures_1-3\\170222_rasmc_active_pol2_0v24_meta.pdf',width = 10,height = 8)
plot(1:320, pol2_0_meta,type='l',col='red',ylim =c(0,2.5),xaxt='n',xlab='',ylab='rpm/bp',main='active genes 0H (red) vs. 24H (black)')
lines(1:320, pol2_24_meta,type='l',col='black')
axis(1,c(0,60,260,320),c('-3kb','TSS','END','+3kb'))
dev.off()

pdf(file='C:\\Users\\rhirsch\\Documents\\rasmc_figures_1-3\\170222_rasmc_active_pol2_2v2J_meta.pdf',width = 10,height = 8)
plot(1:320, pol2_2_meta,type='l',col='red',ylim =c(0,2.5),xaxt='n',xlab='',ylab='rpm/bp',main='active genes 2H (red) vs. 2H+JQ1 (black)')
lines(1:320, pol2_2J_meta,type='l',col='black')
axis(1,c(0,60,260,320),c('-3kb','TSS','END','+3kb'))
dev.off()

pdf(file='C:\\Users\\rhirsch\\Documents\\rasmc_figures_1-3\\170222_rasmc_active_pol2_24v24J_meta.pdf',width = 10,height = 8)
plot(1:320, pol2_24_meta,type='l',col='red',ylim =c(0,5),xaxt='n',xlab='',ylab='rpm/bp',main='active genes 24H (red) vs. 24H+JQ1 (black)')
lines(1:320, pol2_24J_meta,type='l',col='black')
axis(1,c(0,60,260,320),c('-3kb','TSS','END','+3kb'))
dev.off()

###########################################
############Gained/Lost Genes Metas 24H####
###########################################

lost_genes=c()
for(i in 1:1000){
  row = which(as.character(active_genes[,2]) == as.character(orderTable[i,1]))
  lost_genes = c(lost_genes,row)
}

gained_genes = c()
for(i in 7726:8725){
  row = which(as.character(active_genes[,2])==as.character(orderTable[i,1]))
  gained_genes = c(gained_genes,row)
}

pol2_rows_up = c()
for(i in 1:length(gained_genes)){
  row = which(rownames(pol2_0)==active_genes[i,1])
  pol2_rows_up = c(pol2_rows_up,row)
}

pol2_rows_down = c()
for(i in 1:length(lost_genes)){
  row = which(rownames(pol2_0) == active_genes[i,1])
  pol2_rows_down = c(pol2_rows_down,row)
}


pol2_24_meta_up = apply(pol2_24[pol2_rows_up,],2,mean,na.rm=TRUE)
print(pol2_24_meta_up[1:5])
print(pol2_24_meta_up[150:155])
print(pol2_24_meta_up[315:320])
pol2_24_meta_down = apply(pol2_24[pol2_rows_down,],2,mean,na.rm=TRUE)
print(pol2_24_meta_down[1:5])
print(pol2_24_meta_down[150:155])
print(pol2_24_meta_down[315:320])


pdf(file='C:\\Users\\rhirsch\\Documents\\rasmc_figures_1-3\\170222_rasmc_gained_vs_lost_24H_meta_750.pdf',width = 10,height = 8)
plot(1:320, pol2_24_meta_up,type='l',col='red',ylim =c(0,2.5),xaxt='n',xlab='',ylab='rpm/bp',main='lost genes 24H (black) vs. gained genes 24H (red)')
lines(1:320, pol2_24_meta_down,type='l',col='black')
axis(1,c(0,60,260,320),c('-3kb','TSS','END','+3kb'))
dev.off()

###########################################
############Gained/Lost Genes Metas 2H#####
###########################################
table2_genes = as.character(rownames(orderTable2))
lost_genes2=c()
for(i in 1:500){
  row = which(as.character(active_genes[,2]) == table2_genes[i])
  lost_genes2 = c(lost_genes2,row)
}

gained_genes2 = c()
for(i in 7726:8725){
  row = which(as.character(active_genes[,2])== table2_genes[i])
  gained_genes2 = c(gained_genes2,row)
}

pol2_rows_up2 = c()
for(i in 1:length(gained_genes2)){
  row = which(rownames(pol2_0)==active_genes[i,1])
  pol2_rows_up2 = c(pol2_rows_up2,row)
}

pol2_rows_down2 = c()
for(i in 1:length(lost_genes2)){
  row = which(rownames(pol2_0) == active_genes[i,1])
  pol2_rows_down2 = c(pol2_rows_down2,row)
}


pol2_2_meta_up = apply(pol2_2[pol2_rows_up2,],2,mean,na.rm=TRUE)
pol2_2_meta_down = apply(pol2_2[pol2_rows_down2,],2,mean,na.rm=TRUE)



pdf(file='C:\\Users\\rhirsch\\Documents\\rasmc_figures_1-3\\170222_rasmc_gained_vs_lost_2H_meta.pdf',width = 10,height = 8)
plot(1:320, pol2_2_meta_up,type='l',col='red',ylim =c(0,2.5),xaxt='n',xlab='',ylab='rpm/bp',main='lost genes 2H (black) vs. gained genes 2H (red)')
lines(1:320, pol2_2_meta_down,type='l',col='black')
axis(1,c(0,60,260,320),c('-3kb','TSS','END','+3kb'))
dev.off()



