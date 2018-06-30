error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

library(plotrix) #<- needed for standard error std.error function

#load in enhancerPromoter tables

a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_UNSTIM_REP1_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_UNSTIM_REP1_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_2H_REP2_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_2H_REP2_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_24H_REP2_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_24H_REP2_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)


b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_UNSTIM_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_UNSTIM_NEW_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_24H_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_24H_NEW_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)


#create matrix with preserved order

a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE

#create tables for promoter data, enhancer data, and cumulative data for each sample

promoterTable = cbind(a0H_GENE_TABLE[,2],b0H_GENE_TABLE[,2],a2H_GENE_TABLE[,2],b2H_GENE_TABLE[,2],a24H_GENE_TABLE[,2],b24H_GENE_TABLE[,2])
#print(promoterTable[1:10,])
enhancerTable = cbind(a0H_GENE_TABLE[,3],b0H_GENE_TABLE[,3],a2H_GENE_TABLE[,3],b2H_GENE_TABLE[,3],a24H_GENE_TABLE[,3],b24H_GENE_TABLE[,3])
#print(enhancerTable[1:10,])
cumulativeTable = cbind((promoterTable[,1]+enhancerTable[,1]),(promoterTable[,2]+enhancerTable[,2]),(promoterTable[,3]+enhancerTable[,3]),
                        (promoterTable[,4]+enhancerTable[,4]),(promoterTable[,5]+enhancerTable[,5]),(promoterTable[,6]+enhancerTable[,6]))
#print(cumulativeTable[1:10,])

#turn tables into matrices
p_matrix = as.matrix(promoterTable)
e_matrix = as.matrix(enhancerTable)
c_matrix = as.matrix(cumulativeTable)

colnames(p_matrix)=c('UNSTIM_REP1','UNSTIM_NEW','2H_REP2','2H_NEW','24H_REP2','24H_NEW')
colnames(e_matrix)=c('UNSTIM_REP1','UNSTIM_NEW','2H_REP2','2H_NEW','24H_REP2','24H_NEW')
colnames(c_matrix)=c('UNSTIM_REP1','UNSTIM_NEW','2H_REP2','2H_NEW','24H_REP2','24H_NEW')

rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE


#find data for nhr and jf modules
NHR_list = c('Thra','Thrb', 'Nr1d1', 'Rara', 'Rxra', 'Rxrb', 'Rarg', 'Esrra', 'Nr2f6','Nr2f2','Ppara')
JF_list = c('Jun', 'Junb', 'Jund', 'Fos','Fosl1', 'Fosl2')

nhr_rows=c()
jf_rows=c()

for(name in NHR_list){
  row=which(rownames(p_matrix)==name)
  nhr_rows=c(row,nhr_rows)
}

for(name in  JF_list){
  row=which(rownames(p_matrix)==name)
  jf_rows=c(row,jf_rows)
}


p_jf=p_matrix[jf_rows,]
e_jf=e_matrix[jf_rows,]
c_jf=c_matrix[jf_rows,]

p_nhr=p_matrix[nhr_rows,]
e_nhr=e_matrix[nhr_rows,]
c_nhr=c_matrix[nhr_rows,]

c_table=rbind(c_jf,c_nhr)

#median normalize columns (samples) to account for replicate discrepancies
med_unstim_old=median(c_matrix[,1])
med_2_old=median(c_matrix[,3])
med_24_old=median(c_matrix[,5])
med_unstim_new=median(c_matrix[,2])
med_2_new=median(c_matrix[,4])
med_24_new=median(c_matrix[,6])

#create matrix of log2 fold change of cumulative data medians
medianMatrix =cbind(log2(c_table[,1]/med_unstim_old),log2(c_table[,2]/med_unstim_new),log2(c_table[,3]/med_2_old),
                    log2(c_table[,4]/med_2_new),log2(c_table[,5]/med_24_old),log2(c_table[,6]/med_24_new))

colnames(medianMatrix)=c('UNSTIM_REP1','UNSTIM_NEW','2H_REP2','2H_NEW','24H_REP2','24H_NEW')

#take average across log2 FC median normalized replicates

avg_0 = (medianMatrix[,1]+medianMatrix[,2])/2
avg_2 = (medianMatrix[,3]+medianMatrix[,4])/2
avg_24 = (medianMatrix[,5]+medianMatrix[,6])/2


#avg_med = cbind(avg_0,avg_2,avg_24)
#t_avg = t(avg_med)

change_2v0 = avg_2-avg_0
change_24v0 = avg_24-avg_0
change_24v2 = avg_24-avg_2


#fc_matrix = cbind(change_2v0,change_24v0,change_24v2)

#create average fold change matrix for plotting exemplars at 0h vs 0h, 2h vs 0h, and 24h vs 0h
avg_fc_matrix=cbind((avg_0-avg_0),change_2v0,change_24v0)

#t_fc = t(fc_matrix)
#nhr_t_fc=t_fc[,6:15]


#########################################################################
###########################EXEMPLAR PLOTS################################
#########################################################################





moduleMatrix = avg_fc_matrix[1:5,]
colnames(moduleMatrix)=c('0H','2H','24H')

#Plot jf exemplar plots
pdf('/storage/cylin/grail/projects/rasmc_all/figures/jun_fos_median_normalized_exemplar_plots.pdf')

  
plot(1:ncol(moduleMatrix),rep(0,ncol(moduleMatrix[,])),ylim =c(-.5,.5),cex=0,
       xlab='Samples',xaxt='n',ylab='log2 fold change median normalized Brd4 AUC at promoter and enhancer vs 0H',main ='Jun/Fos')
axis(1,1:ncol(moduleMatrix),colnames(moduleMatrix[,]),las=2)
for(i in 1:nrow(moduleMatrix)){
  lines(1:ncol(moduleMatrix),moduleMatrix[i,], col = rgb(0.75,0.75,0.75,0.5),lwd=1)
}
lines(1:ncol(moduleMatrix),apply(moduleMatrix,2,mean),col='black',lwd=4)
error.bar(1:ncol(moduleMatrix),apply(moduleMatrix,2,mean),apply(moduleMatrix,2,std.error)) #standard error of the mean
#error.bar(1:ncol(moduleMatrix),apply(moduleMatrix,2,mean),apply(moduleMatrix,2,sd)) #sd

legend(1,quantile(moduleMatrix,.90),paste("n =",nrow(moduleMatrix)))



dev.off()



####################################################################

moduleMatrix = avg_fc_matrix[6:15,]
colnames(moduleMatrix)=c('0H','2H','24H')

#Plot nhr exemplar plots
pdf('/storage/cylin/grail/projects/rasmc_all/figures/nhr_median_normalized_exemplar_plots.pdf')


plot(1:ncol(moduleMatrix),rep(0,ncol(moduleMatrix[,])),ylim =c(-.5,.5),cex=0,
     xlab='Samples',xaxt='n',ylab='log2 fold change median normalized Brd4 AUC at promoter and enhancer vs 0H',main ='NHR')
axis(1,1:ncol(moduleMatrix),colnames(moduleMatrix[,]),las=2)
for(i in 1:nrow(moduleMatrix)){
  lines(moduleMatrix[i,], col = rgb(1,.5,0,0.5),lwd=1)
}
lines(1:ncol(moduleMatrix),apply(moduleMatrix,2,mean),col='red',lwd=4)
error.bar(1:ncol(moduleMatrix),apply(moduleMatrix,2,mean),apply(moduleMatrix,2,std.error),col='red') #standard error of the mean

legend(1,quantile(moduleMatrix,.90),paste("n =",nrow(moduleMatrix)))



dev.off()






# #========================================================================
# #============================DOT PLOTS===================================
# #========================================================================
# dotPlotGenes <- function(geneList,expressionTable){
#   
#   expMatrix = matrix(nrow = length(geneList),ncol=10)
#   #colorVector = brewer.pal(11,'Spectral')
#   colorVector = colorRampPalette(c('red','grey'))(length(geneList))
#   for(i in 1:length(geneList)){
#     geneName = geneList[i]
#     geneRow = which(row.names(expressionTable) == geneName)
#     if(length(geneRow)==0){print(geneName)}
#     expMatrix[i,] = as.numeric(expressionTable[geneRow,])
#   }
#   
#   dotMatrix = expMatrix
#   
#   dotMatrix[which(dotMatrix< -1)] = -1
#   
#   xRange = c(0,(nrow(expMatrix)))
#   yRange = c(-1,1)
#   
#   plot(0,0,ylim = yRange,xlim = xRange,cex=0,xaxt='n',ylab = 'log2 fold change median normalized Brd4 signal',xlab='',yaxt='n',las=2)
#   axis(1,seq(0.5,nrow(expMatrix)-0.5,1),labels = geneList,las=2)
#   axis(2,c(-1,0,1))
#   for(i in 1:nrow(expMatrix)){
#     
#     x = i - 0.5
#     xVector = jitter(rep(x,ncol(expMatrix)),amount =0.1)
#     points(xVector,dotMatrix[i,],col = colorVector[i],pch=16,cex=1)
#     segments(x-.25,mean(expMatrix[i,]),x+.25,mean(expMatrix[i,]),lwd=4,col = rgb(0.5,0.5,0.5,0.5))
#     
#     
#   }
#   
#   
#   
# }
# 
# names=c('change_2v0','change_24v0','change_24v2')
# 
# pdf(file='C:/Users/rhirsch/Documents/nhr_module_brd4_fc_dot_plots.pdf',width = 10,height =6)
# dotPlotGenes(names,t_fc[,6:15])
# dev.off()
# 
# names=c('avg_0','avg_2', 'avg_24')
# pdf(file='C:/Users/rhirsch/Documents/nhr_module_brd4_med_dot_plots.pdf',width = 10,height =6)
# dotPlotGenes(names,t_avg[,6:15])
# dev.off()




# #========================================================================
# #============================DOT PLOTS===================================
# #========================================================================
# 
# dotPlotGenes <- function(geneList,expressionTable){
#   
#   expMatrix = matrix(nrow = length(geneList),ncol=5)
#   #colorVector = brewer.pal(11,'Spectral')
#   colorVector = colorRampPalette(c('red','grey'))(length(geneList))
#   for(i in 1:length(geneList)){
#     geneName = geneList[i]
#     geneRow = which(row.names(expressionTable) == geneName)
#     if(length(geneRow)==0){print(geneName)}
#     expMatrix[i,] = as.numeric(expressionTable[geneRow,])
#   }
#   
#   dotMatrix = expMatrix
#   
#   dotMatrix[which(dotMatrix< -1)] = -1
#   
#   xRange = c(0,(nrow(expMatrix)))
#   yRange = c(-1,1)
#   
#   plot(0,0,ylim = yRange,xlim = xRange,cex=0,xaxt='n',ylab = 'log2 fold change median normalized Brd4 signal',xlab='',yaxt='n',las=2)
#   axis(1,seq(0.5,nrow(expMatrix)-0.5,1),labels = geneList,las=2)
#   axis(2,c(-1,0,1))
#   for(i in 1:nrow(expMatrix)){
#     
#     x = i - 0.5
#     xVector = jitter(rep(x,ncol(expMatrix)),amount =0.1)
#     points(xVector,dotMatrix[i,],col = colorVector[i],pch=16,cex=1)
#     segments(x-.25,mean(expMatrix[i,]),x+.25,mean(expMatrix[i,]),lwd=4,col = rgb(0.5,0.5,0.5,0.5))
#     
#     
#   }
#   
#   
#   
# }
# 
# names=c('change_2v0','change_24v0','change_24v2')
# 
# pdf(file='C:/Users/rhirsch/Documents/junfos_module_brd4_fc_dot_plots.pdf',width = 10,height =6)
# dotPlotGenes(names,t_fc[,1:5])
# dev.off()
# 
# names=c('avg_0','avg_2', 'avg_24')
# pdf(file='C:/Users/rhirsch/Documents/junfos_module_brd4_med_dot_plots.pdf',width = 10,height =6)
# dotPlotGenes(names,t_avg[,1:5])
# dev.off()
# 
# 
# t.test(t_fc[1,1:5],t_fc[2,1:5])
# t.test(t_fc[1,1:5],t_fc[3,1:5])
# t.test(t_fc[2,1:5],t_fc[3,1:5])
# 
# t.test(t_fc[1,6:15],t_fc[2,6:15])
# t.test(t_fc[1,6:15],t_fc[3,6:15])
# t.test(t_fc[2,6:15],t_fc[3,6:15])
# 
# 
# 
# ######################################################################################################
# ######################################################################################################
# 
# spry2_row = which(rownames(c_matrix)=='Spry2')
# col1a1_row = which(rownames(c_matrix)=='Col1a1')
# 
# sp_col_table= cbind(c_matrix[spry2_row,],c_matrix[col1a1_row,])
# 
# colnames(sp_col_table)=c('Spry2','Col1a1')
# 
# med_spc_table = rbind(log2(sp_col_table[1,]/med_unstim_old),log2(sp_col_table[2,]/med_unstim_new),log2(sp_col_table[3,]/med_2_old),
#                       log2(sp_col_table[4,]/med_2_new),log2(sp_col_table[5,]/med_24_old),log2(sp_col_table[6,]/med_24_new))
# 
# avg_spc_0 = (med_spc_table[1,]+med_spc_table[2,])/2
# avg_spc_2 = (med_spc_table[3,]+med_spc_table[4,])/2
# avg_spc_24 = (med_spc_table[5,]+med_spc_table[2,])/2


