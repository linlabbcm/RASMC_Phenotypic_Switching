#=========================================================
#========================HELPER FUNCTIONS=================
#=========================================================

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


compareVectors <- function(v1,v2,name1,name2,title,nBins = 100,nIter = 1000,yMin='',yMax=''){
  
  #first one as a line graph
  
  v1Order = order(v1)
  colorSpectrum <- colorRampPalette(c("green","black","black","red"))(100)
  
  #setting a color data range
  minValue <- -2
  maxValue <- 2
  color_cuts <- seq(minValue,maxValue,length=100)
  color_cuts <- c(min(v1,na.rm=TRUE), color_cuts,max(v1,na.rm=TRUE))
  
  
  #add one extra min color to even out sampling
  colorSpectrum <- c(colorSpectrum[1],colorSpectrum[1],colorSpectrum)
  
  colorVector = c()
  for(i in v1Order){
    delta = v1[i]
    color = colorSpectrum[max(which(color_cuts <= delta))]
    colorVector =c(colorVector,color)
  }
  
  #set the layout
  m =matrix(data = c(1,1,2,2,2,2),ncol=1)
  layout(m)
  
  
  plot(1:length(v1),v1[v1Order],ylim =c(-1.2*max(abs(v1)),1.2*max(abs(v1))),type='l',lwd=2,ylab=name1,xlab=paste(c('Regions ranked by',name1),sep=' '))
  lines(1:length(v1),v1[v1Order],type='h',lwd=3,col=colorVector)
  abline(h=0)
  
  #now for variable #2
  
  #matrix to bind the data
  binSize = length(v1)/nBins
  binMatrix = matrix(ncol = nBins,nrow=binSize)
  i = 1
  for(i in 1:nBins){
    start = 1 + (i-1)*binSize
    stop = i*binSize
    
    binMatrix[,i] = as.numeric(v2[v1Order[start:stop]])
    
    
  }
  
  meanMatrix = matrix(ncol = nBins,nrow=nIter)
  
  for(i in 1:nIter){
    for(j in 1:nBins){
      meanMatrix[i,j] = mean(sample(binMatrix[,j],binSize,replace=TRUE))
    }
  }
  
  meanVector = apply(meanMatrix,2,mean)
  upperError = apply(meanMatrix,2,quantile,probs=0.975) - apply(meanMatrix,2,mean)
  lowerError = apply(meanMatrix,2,mean) - apply(meanMatrix,2,quantile,probs=0.025)
  
  if(yMax == ''){
    yMax = 1.02* max(apply(meanMatrix,2,quantile,probs=0.975))
  }
  
  if(yMin == ''){
    yMin = 0.98 * min(apply(meanMatrix,2,quantile,probs=0.025))
  }
  plot(1:nBins,meanVector,ylim =c(yMin,yMax),ylab=name2,xlab='',cex=1,col=rgb(0.2,0.2,0.2,0.2),pch=16,xaxt='n',main=title)
  error.bar(1:nBins,meanVector,upperError,lowerError,length =0.01,col=rgb(0.2,0.2,0.2,0.2))	
  
  x=1:nBins
  lw1 = loess(meanVector~x)
  lines(x,lw1$fitted,col='blue',lwd=2)
  
}



#example
# pdf(file = "/home/rahirsch/Desktop/test_pdf.pdf")
# v1 = jitter(1:1000,amount = 100) #makes a jittery vector that goes from 1 to 1000
# v2 = jitter((1:1000)^2,amount = 100) #makes a jittery exponential vector
# 
# compareVectors(v1,v2,'v1','v2','example')
# dev.off()

a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_UNSTIM_REP1_0_STITCH_-_JQ1/RASMC_H3K27AC_UNSTIM_REP1_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_2H_REP2_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_2H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_24H_REP2_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_24H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)

b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_UNSTIM_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_UNSTIM_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_2H_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_2H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_24H_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_24H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)

a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE

enhancerTableH = cbind((a0H_GENE_TABLE[,2]+b0H_GENE_TABLE[,2])/2,(a2H_GENE_TABLE[,2]+b2H_GENE_TABLE[,2])/2,(a24H_GENE_TABLE[,2]+b24H_GENE_TABLE[,2])/2)
promoterTableH = cbind((a0H_GENE_TABLE[,3]+b0H_GENE_TABLE[,3])/2,(a2H_GENE_TABLE[,3]+b2H_GENE_TABLE[,3])/2,(a24H_GENE_TABLE[,3]+b24H_GENE_TABLE[,3])/2)
cumulativeTableH = cbind((enhancerTableH[,1]+promoterTableH[,1]),(enhancerTableH[,2]+promoterTableH[,2]),(enhancerTableH[,3]+promoterTableH[,3]))

c0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_UNSTIM_REP1_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_UNSTIM_REP1_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
c2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_2H_REP2_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_2H_REP2_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
c24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_24H_REP2_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_24H_REP2_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
d0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_UNSTIM_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_UNSTIM_NEW_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
d2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)
d24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_24H_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_24H_NEW_0_STITCH_-_JQ1_h3k_gff_GENE_TABLE.txt', sep = '\t', header = TRUE)

promoterTableB = cbind((c0H_GENE_TABLE[,2]+d0H_GENE_TABLE[,2])/2,(c2H_GENE_TABLE[,2]+d2H_GENE_TABLE[,2])/2,(c24H_GENE_TABLE[,2]+d24H_GENE_TABLE[,2])/2)
#print(promoterTable[1:10,])
enhancerTableB = cbind((c0H_GENE_TABLE[,3]+d0H_GENE_TABLE[,3])/2,(c2H_GENE_TABLE[,3]+d2H_GENE_TABLE[,3])/2,(c24H_GENE_TABLE[,3]+d24H_GENE_TABLE[,3])/2)
#print(enhancerTable[1:10,])
cumulativeTableB = cbind((promoterTableB[,1]+enhancerTableB[,1]),(promoterTableB[,2]+enhancerTableB[,2]),(promoterTableB[,3]+enhancerTableB[,3]))
#print(cumulativeTable[1:10,])


enhPro_table = cbind(cumulativeTableH[,1:3],cumulativeTableB[,1:3])

colnames(enhPro_table) = c('H3K27AC_UNSTIM', 'H3K27AC_2H', 'H3K27AC_24H', 'BRD4_UNSTIM','BRD4_2H', 'BRD4_24H')
rownames(enhPro_table) = a0H_GENE_TABLE[,1]

rna_table = read.delim('/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_exprs_fpkm_means.txt', header=TRUE, sep = '\t')

genes=rownames(enhPro_table)

####check to see if all genes and rna overlap###


# geneRows = c()
# for(i in 1:length(genes)){
#   geneRow = which(rownames(rna_table)==genes[i])
#   geneRows = c(geneRows,geneRow)
#   }

zeroH0h = which(enhPro_table[,1] != 0)
zeroH2h = which(enhPro_table[,2] !=0)
zeroH24h = which(enhPro_table[,3] != 0)
zeroB0h = which(enhPro_table[,4] != 0)
zeroB2h = which(enhPro_table[,5] !=0)
zeroB24h = which(enhPro_table[,6] != 0)

zeroFree = intersect(zeroH0h,zeroH24h)
zeroFree = intersect(zeroFree,zeroB0h)
zeroFree = intersect(zeroFree,zeroB24h)
zeroFree = intersect(zeroFree,zeroB2h)
zeroFree = intersect(zeroFree,zeroH2h)

cropTable = enhPro_table[zeroFree,]
cropGenes = rownames(cropTable)

exp_rows = c()
for(i in 1:length(cropGenes)){
  exp_row = which(rownames(rna_table) == cropGenes[i])
  exp_rows = c(exp_rows,exp_row)
}

cropTable = cbind(cropTable, rna_table[exp_rows,1],rna_table[exp_rows,2],rna_table[exp_rows,4])

colnames(cropTable) = c('H3K27AC_UNSTIM','H3K27AC_2H', 'H3K27AC_24H', 'BRD4_UNSTIM', 'BRD4_2H','BRD4_24H', 'RNA_0H', 'RNA_2H','RNA_24H')
#############ADD LOG COLUMNS TO SCRIPTS

logMatrix = log2(cropTable)


#############MAKE COMP VECTORS FIGURES
# 
#double fold change BRD4 24 vs H3k 24  hour
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/double_fold_change_brd4_24_vs_h3k_24_9037.pdf")
v1 = logMatrix[,6] - logMatrix[,4]
v2 = logMatrix[,3] - logMatrix[,1]

compareVectors(v1,v2,'Log2 Fold Change 24H BRD4','Log2 Fold Change 24H H3k27ac', '')
dev.off()

#double fold change BRD4 24H vs exp 24 hour

pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/double_fold_change_brd4_24_vs_rna_24_9037.pdf")
v1 = logMatrix[,6] - logMatrix[,4]
v2 = logMatrix[,9] - logMatrix[,7]

compareVectors(v1,v2,'Log2 Fold Change 24H BRD4','Log2 Fold Change 24H Expression', '')
dev.off()

#double fold change BRD4 2 vs H3k 2  hour
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/double_fold_change_brd4_2_vs_h3k_2_9037.pdf")
v1 = logMatrix[,5] - logMatrix[,4]
v2 = logMatrix[,2] - logMatrix[,1]

compareVectors(v1,v2,'Log2 Fold Change 2H BRD4','Log2 Fold Change 2H H3k27ac', '')
dev.off()

#double fold change BRD4 2 vs exp 2 hour

pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/double_fold_change_brd4_2_vs_rna_2_9037.pdf")
v1 = logMatrix[,5] - logMatrix[,4]
v2 = logMatrix[,8] - logMatrix[,7]

compareVectors(v1,v2,'Log2 Fold Change 2H BRD4','Log2 Fold Change 2H Expression', '')
dev.off()

# 
# # #double fold change 24 vs 2 hour rna
# pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/double_fold_change_h3k27ac_0h_24h_vs_2h_rna_log2_change.pdf")
# v1 = logMatrix[,6] - logMatrix[,4]
# v2 = logMatrix[,2] - logMatrix[,1]
# 
# compareVectors(v1,v2,'Log2 Fold Change 24H Brd4','Log2 Fold Change 2H Expression', '')
# dev.off()
# 
# #double fold change 24 vs 24 hour
# 
# pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/double_fold_change_h3k27ac_0h_24h_vs_24h_rna_log2_change.pdf")
# v1 = logMatrix[,6] - logMatrix[,4]
# v2 = logMatrix[,3] - logMatrix[,1]
# 
# compareVectors(v1,v2,'Log2 Fold Change 24H Brd4','Log2 Fold Change 24H Expression', '')
# dev.off()

# #2 vs 24 hour brd4 h3k27ac
#
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/compare_vectors_brd4_2h_vs_24h_h3k27ac_log2_change.pdf")
v1 = logMatrix[,6] - logMatrix[,5]
v2 = logMatrix[,3] - logMatrix[,2]
#
compareVectors(v1,v2,'Log2 Fold Change 2v24H Brd4','Log2 Fold Change 2v24H H3k27ac', '')
dev.off()

# #2 vs 24 hour brd4 rna
#
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/compare_vectors_brd4_2h_vs_24h_rna_log2_change.pdf")
v1 = logMatrix[,6] - logMatrix[,5]
v2 = logMatrix[,9] - logMatrix[,8]
#
compareVectors(v1,v2,'Log2 Fold Change 2v24H Brd4','Log2 Fold Change 2v24H Expression', '')
dev.off()











#######MAKING VECTOR TABLE##############
doubleLogMatrixH3k = logMatrix - logMatrix[,1]
doubleLogMatrixBRD4 = logMatrix-logMatrix[,4]
doubleLogMatrixRNA = logMatrix - logMatrix[,7]

foo = order(doubleLogMatrixBRD4[,6])

names = rownames(doubleLogMatrixBRD4[foo,])

orderMatrix = cbind(cropTable[foo,4],cropTable[foo,5], cropTable[foo,6], cropTable[foo,1],cropTable[foo,2], cropTable[foo,3],cropTable[foo,7], cropTable[foo,8],cropTable[foo,9],doubleLogMatrixBRD4[foo,5],doubleLogMatrixBRD4[foo,6],doubleLogMatrixH3k[foo,2], doubleLogMatrixH3k[foo,3], doubleLogMatrixRNA[foo,8], doubleLogMatrixRNA[foo,9])
columns = c('brd4_unstim', 'brd4_2', 'brd4_24','h3k27ac_unstim','h3k27ac_2h','h3k27ac_24h','rna_unstim','rna_2h','rna_24h', 'brd4_log2_2_v_0', 'brd4_log2_24_v_0','h3k27ac_log2_2_v_0', 'h3k27ac_log2_24_v_0', 'rna_log2_2_v_0','rna_log2_24_v_0')
colnames(orderMatrix) = columns
rownames(orderMatrix) = names


derpMatrix = orderMatrix

#write.table(orderMatrix,file='/storage/cylin/grail/projects/rasmc_all/TF_Network/170206_waterfall_genes_log2_24_v_0_ordered_BRD4_H3k_RNA_0_2_24_and_log2.txt',quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)



doubleLogMatrixH3k = logMatrix - logMatrix[,1]
doubleLogMatrixBRD4 = logMatrix-logMatrix[,4]
doubleLogMatrixRNA = logMatrix - logMatrix[,7]

foop = order(doubleLogMatrixBRD4[,5])

names = rownames(doubleLogMatrixBRD4[foop,])

orderMatrix2 = cbind(cropTable[foop,4],cropTable[foop,5], cropTable[foop,6], cropTable[foop,1],cropTable[foop,2], cropTable[foop,3],cropTable[foop,7], cropTable[foop,8],cropTable[foop,9],doubleLogMatrixBRD4[foop,5],doubleLogMatrixBRD4[foop,6],doubleLogMatrixH3k[foop,2], doubleLogMatrixH3k[foop,3], doubleLogMatrixRNA[foop,8], doubleLogMatrixRNA[foop,9])
columns = c('brd4_unstim', 'brd4_2', 'brd4_24','h3k27ac_unstim','h3k27ac_2h','h3k27ac_24h','rna_unstim','rna_2h','rna_24h', 'brd4_log2_2_v_0', 'brd4_log2_24_v_0','h3k27ac_log2_2_v_0', 'h3k27ac_log2_24_v_0', 'rna_log2_2_v_0','rna_log2_24_v_0')
colnames(orderMatrix2) = columns
rownames(orderMatrix2) = names


derpMatrix = orderMatrix2

#write.table(orderMatrix2,file='/storage/cylin/grail/projects/rasmc_all/TF_Network/170206_waterfall_genes_log2_2_v_0_ordered_BRD4_H3k_RNA_0_2_24_and_log2.txt',quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)


doubleLogMatrixH3k = logMatrix - logMatrix[,2]
doubleLogMatrixBRD4 = logMatrix-logMatrix[,5]
doubleLogMatrixRNA = logMatrix - logMatrix[,8]


foop = order(doubleLogMatrixBRD4[,6])

names = rownames(doubleLogMatrixBRD4[foop,])

orderMatrix2v24 = cbind(cropTable[foop,4],cropTable[foop,5], cropTable[foop,6], cropTable[foop,1],cropTable[foop,2], cropTable[foop,3],cropTable[foop,7], cropTable[foop,8],cropTable[foop,9],doubleLogMatrixBRD4[foop,6],doubleLogMatrixH3k[foop,3],doubleLogMatrixRNA[foop,9])
columns = c('brd4_unstim', 'brd4_2', 'brd4_24','h3k27ac_unstim','h3k27ac_2h','h3k27ac_24h','rna_unstim','rna_2h','rna_24h', 'brd4_log2_2_v_24','h3k27ac_log2_2_v_24', 'rna_log2_2_v_24')
colnames(orderMatrix2v24) = columns
rownames(orderMatrix2v24) = names


derpMatrix = orderMatrix2v24

write.table(orderMatrix2v24,file='/storage/cylin/grail/projects/rasmc_all/TF_Network/170206_waterfall_genes_log2_2_v_24_ordered_BRD4_H3k_RNA_0_2_24_and_log2.txt',quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)

