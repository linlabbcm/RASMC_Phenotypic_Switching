expTable = read.delim('/storage/cylin/grail/projects/rasmc/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_all_fpkm_exprs_norm.txt')

diffGenes = read.delim('/storage/cylin/grail/projects/rasmc_all/STRING/known_interaction_proteins.txt')

diffGenesList = as.character(diffGenes[,1])

diffGenesList = sort(diffGenesList)

sillyColumns = setdiff(1:ncol(expTable),grep('JQ1',colnames(expTable)))
expMatrix = as.matrix(expTable[,sillyColumns])


diffMatrix = matrix(nrow=length(diffGenesList),ncol=ncol(expMatrix))

rownames(diffMatrix) = diffGenesList
colnames(diffMatrix) = colnames(expMatrix)


for(i in 1:length(diffGenesList)){
  
  geneName = diffGenesList[i]
  expRow = which(rownames(expMatrix)==geneName)
  diffMatrix[i,] = expMatrix[expRow,]
  
}


write.table(diffMatrix, file = "/storage/cylin/grail/projects/rasmc_all/tables/known_int_heatmap_matrix.txt", append = FALSE, 
            quote = FALSE, sep = "\t", eol = '\n', na = "NA", dec = '.', row.names = FALSE, col.names = TRUE,
            qmethod = c("escape", "double"), fileEncoding = "")


#now we want to make a log2 row median normalized expression matrix

medianVector = apply(diffMatrix,1,median) 

medianMatrix = log2(diffMatrix/medianVector)

#zero hour mean normalization

#meanVector = apply(diffMatrix[,1:2],1,mean)

#meanMatrix = log2(diffMatrix)-log2(meanVector)

#meanMatrix <- meanMatrix[,-1]

#===================================================================
#======================CLUSTERING EXPRESSION========================
#===================================================================

expCorDist = as.dist(1-cor(t(medianMatrix)))
expHClust = hclust(expCorDist)

expOrder = expHClust$order
#===================================================================
#=========================MAKING HEATMAPS===========================
#===================================================================


#Set the color spectrum
colorSpectrum <- colorRampPalette(c("blue","white","red"))(100)

#setting a color data range
#minValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.025,names=FALSE)
#maxValue <- quantile(clusterMatrix,na.rm=TRUE,prob=0.975,names=FALSE)

minValue= -4
maxValue = 4


color_cuts <- seq(minValue,maxValue,length=100)
#color_cuts <- seq(-1,1,length=100)
color_cuts <- c(min(medianMatrix), color_cuts,max(medianMatrix)) # this catches the full dynamic range

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum) #this is something stupid in R that you always have to do

#Making png

#clusterPNGFile = paste(outputFolder,genome,'_',analysisName,'_2d_cluster.png',sep='')
png(filename = '/storage/cylin/grail/projects/rasmc_all/figures/170913_expressionHeatmap_known_interactions.png',width = 800,height =800)
layout(matrix(data=c(1,1,1,1,1,2,2),ncol= 7))


image(1:ncol(medianMatrix),1:nrow(medianMatrix),t(medianMatrix[expOrder,]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')


image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Log2 fold vs. median")
dev.off()


#===================================================================
#===========================DENDROGRAMS=============================
#===================================================================

clusterCut <- cutree(expHClust, 10)

##uncomment to plot dendrogram
pdf('/storage/cylin/grail/projects/rasmc_all/figures/known_interactions_dedrogram.pdf')
plot(expHClust, h = -.5)
dev.off()

clusterList <-c(clusterCut)

write.table(clusterList, file = "/storage/cylin/grail/projects/rasmc_all/tables/clusterList_from_known_interactions_heatmap_clusterCut.txt", append = FALSE, 
            quote = FALSE, sep = "\t", eol = '\n', na = "NA", dec = '.', row.names = FALSE, col.names = TRUE,
            qmethod = c("escape", "double"), fileEncoding = "")