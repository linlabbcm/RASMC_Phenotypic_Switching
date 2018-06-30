expTable = read.delim('/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_all_fpkm_exprs_norm.txt')

diffGenes = read.delim('/storage/cylin/grail/projects/rasmc_all/tables/clustergram_genes_1.5_fold_change.txt')

cluster1genes = read.delim('/storage/cylin/grail/projects/rasmc_all/tables/cluster_1_genes_list.txt')
cluster4genes = read.delim('/storage/cylin/grail/projects/rasmc_all/tables/cluster_4_genes_list.txt')
cluster5genes = read.delim('/storage/cylin/grail/projects/rasmc_all/tables/cluster_5_genes_list.txt')
cluster6genes = read.delim('/storage/cylin/grail/projects/rasmc_all/tables/cluster_6_genes_list.txt')

diffGenesList = as.character(diffGenes[,1])

diffGenesList = sort(diffGenesList)

sillyColumns = setdiff(1:ncol(expTable),grep('JQ1',colnames(expTable)))
expMatrix = as.matrix(expTable[,sillyColumns])

expTableJQ1 = cbind(expTable[,1:2],expTable[,6:8], expTable[,12:14])
expMatrixJQ1 = as.matrix(expTableJQ1)

diffMatrix = matrix(nrow=length(diffGenesList),ncol=ncol(expMatrix))
diffMatrixJQ1 = matrix(nrow=length(diffGenesList), ncol=ncol(expMatrixJQ1))


rownames(diffMatrix) = diffGenesList
colnames(diffMatrix) = colnames(expMatrix)


rownames(diffMatrixJQ1) = diffGenesList
colnames(diffMatrixJQ1) = colnames(expMatrixJQ1)


for(i in 1:length(diffGenesList)){
  
  geneName = diffGenesList[i]
  expRow = which(rownames(expMatrix)==geneName)
  diffMatrix[i,] = expMatrix[expRow,]
  
}

for(i in 1:length(diffGenesList)){
  
  geneName = diffGenesList[i]
  expRow = which(rownames(expMatrixJQ1)==geneName)
  diffMatrixJQ1[i,] = expMatrixJQ1[expRow,]
  
}


diffMatrix = cbind(diffMatrix, clusterList)


mean_0H = apply(diffMatrix[,1:2],1,mean)
mean_2H = apply(diffMatrix[,3:5],1,mean)
mean_24H = apply(diffMatrix[,6:8],1,mean)

mean_2Hjq1 = apply(diffMatrixJQ1[,3:5],1,mean)
mean_24Hjq1 = apply(diffMatrixJQ1[,6:8],1,mean)


meanTable = cbind(mean_0H,mean_2H,mean_2Hjq1,mean_24H,mean_24Hjq1)

rownames(meanTable)==rownames(diffMatrix)

clusterArows = c()
clusterBrows = c()
clusterCrows = c()
clusterDrows = c()

for(i in 1:540){
  clusterArow = which(rownames(diffMatrix) == clusterAgenes[i,1])
  clusterArows = append(clusterArows,clusterArow)
}
for(i in 1:712){
  clusterBrow = which(rownames(diffMatrix) == clusterBgenes[i,1])
  clusterBrows = append(clusterBrows,clusterBrow)
}
for(i in 1:1242){
  clusterCrow = which(rownames(diffMatrix) == clusterCgenes[i,1])
  clusterCrows = append(clusterCrows,clusterCrow)
}
for(i in 1:length(clusterDgenes[,1])){
  clusterDrow = which(rownames(diffMatrix) == clusterDgenes[i,1])
  clusterDrows = append(clusterDrows,clusterDrow)
}




for(i in 2:5){meanTable[,i] = log2(meanTable[,i]/meanTable[,1])}

meanTable4 = meanTable[clusterCrows,]

meanTable5 = meanTable[clusterArows,]

meanTable1 = meanTable[clusterBrows,]

meanTable6 = meanTable[clusterDrows,]

pdf('/storage/cylin/grail/projects/rasmc_all/figures/cluster_A_B_C_fold_change_boxplots.pdf')
par(mfrow=c(2,2))

boxplot(meanTable1[,2:5],cex = 0, ylim = c(-2,3), main = 'Cluster B Log2 Fold Change', col = c('white', 'red'))

boxplot(meanTable4[,2:5],cex = 0, ylim = c(-3,2), main = 'Cluster C Log2 Fold Change', col = c('white', 'red'))

boxplot(meanTable5[,2:5],cex = 0, ylim = c(-2,5), main = 'Cluster A Log2 Fold Change', col = c('white', 'red'))

boxplot(meanTable6[,2:5],cex = 0, ylim = c(-2,5), main = 'Cluster D Log2 Fold Change', col = c('white', 'red'))


dev.off()


t.test(meanTable1[,2],meanTable1[,3])
t.test(meanTable1[,4],meanTable1[,5])
t.test(meanTable4[,2],meanTable4[,3])
t.test(meanTable4[,4],meanTable4[,5])
t.test(meanTable5[,2],meanTable5[,3])
t.test(meanTable5[,4],meanTable5[,5])
t.test(meanTable6[,2],meanTable6[,3])
t.test(meanTable6[,4],meanTable6[,5])

