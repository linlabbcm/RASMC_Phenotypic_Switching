a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_UNSTIM_REP1_0_STITCH_-_JQ1/RASMC_H3K27AC_UNSTIM_REP1_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_2H_REP2_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_2H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_JQ1_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_JQ1_2H_REP2_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_JQ1_2H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_24H_REP2_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_24H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_JQ1_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_JQ1_24H_REP2_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_JQ1_24H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)

b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_UNSTIM_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_UNSTIM_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_2H_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_2H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_JQ1_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_JQ1_2H_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_JQ1_2H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_24H_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_24H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_JQ1_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/H3K27AC/RASMC_H3K27AC_PDGF_JQ1_24H_NEW_0_STITCH_-_JQ1/RASMC_H3K27AC_PDGF_JQ1_24H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)






a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE

enhancerTable = cbind((a0H_GENE_TABLE[,2]+b0H_GENE_TABLE[,2])/2,(a2H_GENE_TABLE[,2]+b2H_GENE_TABLE[,2])/2,(a2H_JQ1_GENE_TABLE[,2]+b2H_JQ1_GENE_TABLE[,2])/2,(a24H_GENE_TABLE[,2]+b24H_GENE_TABLE[,2])/2,(a24H_JQ1_GENE_TABLE[,2]+b24H_JQ1_GENE_TABLE[,2])/2)
promoterTable = cbind((a0H_GENE_TABLE[,3]+b0H_GENE_TABLE[,3])/2,(a2H_GENE_TABLE[,3]+b2H_GENE_TABLE[,3])/2,(a2H_JQ1_GENE_TABLE[,3]+b2H_JQ1_GENE_TABLE[,3])/2,(a24H_GENE_TABLE[,3]+b24H_GENE_TABLE[,3])/2,(a24H_JQ1_GENE_TABLE[,3]+b24H_JQ1_GENE_TABLE[,3])/2)
cumulativeTable = cbind((enhancerTable[,1]+promoterTable[,1]),(enhancerTable[,2]+promoterTable[,2]),(enhancerTable[,3]+promoterTable[,3]),(enhancerTable[,4]+promoterTable[,4]),(enhancerTable[,5]+promoterTable[,5]))

e_matrix = as.matrix(enhancerTable)
p_matrix = as.matrix(promoterTable)
c_matrix = as.matrix(cumulativeTable)

rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE

clusterCgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_4_genes_list.txt",header = FALSE)

clusterAgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_5_genes_list.txt",header = FALSE)

clusterBgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_1_genes_list.txt",header = FALSE)

clusterArows = c()
clusterBrows = c()
clusterCrows = c()

for(i in 1:540){
  clusterArow = which(rownames(e_matrix) == clusterAgenes[i,1])
  clusterArows = append(clusterArows,clusterArow)
}
for(i in 1:712){
  clusterBrow = which(a0h_matrix[,1] == clusterBgenes[i,1])
  clusterBrows = append(clusterBrows,clusterBrow)
}
for(i in 1:1242){
  clusterCrow = which(a0h_matrix[,1] == clusterCgenes[i,1])
  clusterCrows = append(clusterCrows,clusterCrow)
}


# table_A <- promoterTable[clusterArows,]
# 
# 
# boxplot(table_A[,1], table_A[,2], table_A[,3],table_A[,4],table_A[,5],ylim=c(0,45000),cex=0, names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
#         ylab= c("Average Signal"), main = c("Promoter Activity","Cluster A"))




pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/CLUSTER_A_B_C_promoter_enhancer_cumulative_boxplots.pdf')
par(mfrow=c(3,1))
#########################################################################################################
table_EA <- e_matrix[clusterArows,]
boxplot(table_EA[,1], table_EA[,2], table_EA[,3],table_EA[,4],table_EA[,5], cex=0, ylim=c(0,20000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Promoter','Cluster A'))

table_EB <- e_matrix[clusterBrows,]
boxplot(table_EB[,1], table_EB[,2], table_EB[,3],table_EB[,4],table_EB[,5], cex=0, ylim=c(0,20000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Promoter','Cluster B'))

table_EC <- e_matrix[clusterCrows,]
boxplot(table_EC[,1], table_EC[,2], table_EC[,3],table_EC[,4],table_EC[,5], cex=0, ylim=c(0,20000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Promoter','Cluster C'))

###########################################################################################################
table_PA <- p_matrix[clusterArows,]
boxplot(table_PA[,1], table_PA[,2], table_PA[,3],table_PA[,4],table_PA[,5], cex=0, ylim=c(0,65000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Enhancer','Cluster A'))

table_PB <- p_matrix[clusterBrows,]
boxplot(table_PB[,1], table_PB[,2], table_PB[,3],table_PB[,4],table_PB[,5], cex=0, ylim=c(0,65000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Enhancer','Cluster B'))

table_PC <- p_matrix[clusterCrows,]
boxplot(table_PC[,1], table_PC[,2], table_PC[,3],table_PC[,4],table_PC[,5], cex=0, ylim=c(0,65000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Enhancer','Cluster C'))

###########################################################################################################3
table_CA <- c_matrix[clusterArows,]
boxplot(table_CA[,1], table_CA[,2], table_CA[,3],table_CA[,4],table_CA[,5], cex=0, ylim=c(0,75000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Cumulative','Cluster A'))

table_CB <- c_matrix[clusterBrows,]
boxplot(table_CB[,1], table_CB[,2], table_CB[,3],table_CB[,4],table_CB[,5], cex=0, ylim=c(0,75000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Cumulative','Cluster B'))

table_CC <- c_matrix[clusterCrows,]
boxplot(table_CC[,1], table_CC[,2], table_CC[,3],table_CC[,4],table_CC[,5], cex=0, ylim=c(0,75000), names = c('0H', '2H','2H+JQ1', '24H','24H+JQ1'), 
        ylab= c("Average Signal"), main = c('Cumulative','Cluster C'))



dev.off()
  


t.test(table_CA[,1], table_CA[,2])
t.test(table_CA[,2], table_CA[,3])
t.test(table_CA[,3], table_CA[,4])
t.test(table_CA[,4], table_CA[,5])

t.test(table_CB[,1], table_CB[,2])
t.test(table_CB[,2], table_CB[,3])
t.test(table_CB[,3], table_CB[,4])
t.test(table_CB[,4], table_CB[,5])

t.test(table_CC[,1], table_CC[,2])
t.test(table_CC[,2], table_CC[,3])
t.test(table_CC[,3], table_CC[,4])
t.test(table_CC[,4], table_CC[,5])

# plotNormGenesSets <-function(A,B,C,cluster_name,output_file){
#   
#   pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/clusters_A_B_C_0h_normalized_enhancer_promoter_cumulative_boxplots.pdf')
#   par(mfrow=c(3,1))
#   # 
#   # for(i in 1:483){
#   #   if(table_EA[i,1] == 0){
#   #     table_EA[i,1] = 1
#   #   }
#   # }
#   
#   #table_EA <- e_matrix[clusterArows,]
#   boxplot(log2(table_EA[,2]/table_EA[,1]), log2(table_EA[,3]/table_EA[,1]),log2(table_EA[,4]/table_EA[,1]),log2(table_EA[,5]/table_EA[,1]), cex=0, ylim=c(-5,10), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster A'))
#   table_EB <- e_matrix[clusterBrows,]
#   boxplot(log2(table_EB[,2]/table_EB[,1]), log2(table_EB[,3]/table_EB[,1]),log2(table_EB[,4]/table_EB[,1]),log2(table_EB[,5]/table_EB[,1]), cex=0, ylim=c(-5,10), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster B'))
#   table_EC <- e_matrix[clusterCrows,]
#   boxplot(log2(table_EC[,2]/table_EC[,1]), log2(table_EC[,3]/table_EC[,1]),log2(table_EC[,4]/table_EC[,1]),log2(table_EC[,5]/table_EC[,1]), cex=0, ylim=c(-5,15), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster C'))
#   
#  #############################################################################################################################################
#   
#   table_PA <- p_matrix[clusterArows,]
#   boxplot(log2(table_PA[,2]/table_PA[,1]), log2(table_PA[,3]/table_PA[,1]),log2(table_PA[,4]/table_PA[,1]),log2(table_PA[,5]/table_PA[,1]), cex=0, ylim=c(-10,10), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster A'))
#   table_PB <- p_matrix[clusterBrows,]
#   boxplot(log2(table_PB[,2]/table_PB[,1]), log2(table_PB[,3]/table_PB[,1]),log2(table_PB[,4]/table_PB[,1]),log2(table_PB[,5]/table_PB[,1]), cex=0,ylim=c(-10,10), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster B'))
#   table_PC <- p_matrix[clusterCrows,]
#   boxplot(log2(table_PC[,2]/table_PC[,1]), log2(table_PC[,3]/table_PC[,1]),log2(table_PC[,4]/table_PC[,1]),log2(table_PC[,5]/table_PC[,1]), cex=0, ylim=c(-10,10), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster C'))
#   
#   ############################################################################################################################################
#   
#   table_EA <- c_matrix[clusterArows,]
#   boxplot(log2(table_CA[,2]/table_CA[,1]), log2(table_CA[,3]/table_CA[,1]),log2(table_CA[,4]/table_CA[,1]),log2(table_CA[,5]/table_CA[,1]), cex=0, ylim=c(-3,7), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster A'))
#   table_EB <- c_matrix[clusterBrows,]
#   boxplot(log2(table_CB[,2]/table_CB[,1]), log2(table_CB[,3]/table_CB[,1]),log2(table_CB[,4]/table_CB[,1]),log2(table_CB[,5]/table_CB[,1]), cex=0, ylim=c(-3,7), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster B'))
#   table_CC <- c_matrix[clusterCrows,]
#   boxplot(log2(table_CC[,2]/table_CC[,1]), log2(table_CC[,3]/table_CC[,1]),log2(table_CC[,4]/table_CC[,1]),log2(table_CC[,5]/table_CC[,1]), cex=0, ylim=c(-3,7), names = c( '2H','2H+JQ1', '24H','24H+JQ1'), 
#           ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster C'))
#   
#   dev.off()
#   
# }
# 
# 
# 
# 
# ###################################################################
# ############# 2H JQ1 vs 2h and 24hjq1 vs 24h#######################
# ###################################################################
# 
# 
# pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/clusters_A_B_C_2h_JQ1_and_24h_JQ1_normalized.pdf')
# par(mfrow=c(3,1))
# 
# #table_EA <- e_matrix[clusterArows,]
# boxplot(log2(table_EA[,3]/table_EA[,2]),log2(table_EA[,5]/table_EA[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Promoter','Cluster A'))
# table_EB <- e_matrix[clusterBrows,]
# boxplot(log2(table_EB[,3]/table_EB[,2]),log2(table_EB[,5]/table_EB[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Promoter','Cluster B'))
# table_EC <- e_matrix[clusterCrows,]
# boxplot(log2(table_EC[,3]/table_EC[,2]),log2(table_EC[,5]/table_EC[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Promoter','Cluster C'))
# 
# #############################################################################################################################################
# 
# table_PA <- p_matrix[clusterArows,]
# boxplot(log2(table_PA[,3]/table_PA[,2]),log2(table_PA[,5]/table_PA[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1','24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Enhancer','Cluster A'))
# table_PB <- p_matrix[clusterBrows,]
# boxplot(log2(table_PB[,3]/table_PB[,2]),log2(table_PB[,5]/table_PB[,4]), cex=0,ylim=c(-3,3), names = c( '2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Enhancer','Cluster B'))
# table_PC <- p_matrix[clusterCrows,]
# boxplot( log2(table_PC[,3]/table_PC[,2]),log2(table_PC[,5]/table_PC[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Enhancer','Cluster C'))
# 
# ############################################################################################################################################
# 
# table_EA <- c_matrix[clusterArows,]
# boxplot(log2(table_CA[,3]/table_CA[,2]),log2(table_CA[,5]/table_CA[,4]), cex=0, ylim=c(-3,3), names = c('2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Cumulative','Cluster A'))
# table_EB <- c_matrix[clusterBrows,]
# boxplot(log2(table_CB[,3]/table_CB[,2]),log2(table_CB[,5]/table_CB[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1','24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Cumulative','Cluster B'))
# table_CC <- c_matrix[clusterCrows,]
# boxplot(log2(table_CC[,3]/table_CC[,2]),log2(table_CC[,5]/table_CC[,4]), cex=0, ylim=c(-3,3), names = c( '2H+JQ1', '24H+JQ1'), 
#         ylab= c("Log2 Fold 2H and 24H Normalized Signals"), main = c('Cumulative','Cluster C'))
# 
# 
# dev.off()
# 

# 
# 
# 
# plotBargraphs <- function(rows, cluster, output_file){
# 
# rows = clusterArows
# cluster = "Cluster A"
# output_file = '/storage/cylin/grail/projects/rasmc_all/figures/cluster_A_enhancer_promoter_cumulative_bargraphs.pdf'  
# 
#   table_E <- e_matrix[rows,]
#   table_P <- p_matrix[rows,]
#   table_C <- c_matrix[rows,]
#   
#   
#   pdf(file = output_file)
#   par(mfrow=c(3,1))
#   for(gene in rownames(table_E)){
#     gene_row = which(rownames(table_E) == gene)
#     centers <- barplot(height = c(table_E[gene_row,1],table_E[gene_row,2],table_E[gene_row,3], table_E[gene_row,4], table_E[gene_row,5]),
#                      beside = true, las = 2,
#                      ylim = c(0,max(table_E[gene_row,1],table_E[gene_row,2],table_E[gene_row,3], table_E[gene_row,4], table_E[gene_row,5])),
#                      cex.names = 0.75, names=c('0H', '2H','2H+JQ1','24H', '24H+JQ1'),
#                      main = c(gene,cluster,"Promoter"),
#                      ylab = "Average Signal",
#                      border = "black", axes = TRUE)
#   
#     gene_row = which(rownames(table_P) == gene)
#     centers <- barplot(height = c(table_P[gene_row,1],table_P[gene_row,2],table_P[gene_row,3], table_P[gene_row,4], table_P[gene_row,5]),
#                        beside = true, las = 2,
#                        ylim = c(0,max(table_P[gene_row,1],table_P[gene_row,2],table_P[gene_row,3], table_P[gene_row,4], table_P[gene_row,5])),
#                        cex.names = .75, names=c('0H', '2H','2H+JQ1','24H', '24H+JQ1'),
#                        main = c(gene,cluster,"Enhancer"),
#                        ylab = "Average Signal",
#                        border = "black", axes = TRUE)
#     
# 
#     gene_row = which(rownames(table_C) == gene)
#     centers <- barplot(height = c(table_C[gene_row,1],table_C[gene_row,2],table_C[gene_row,3], table_C[gene_row,4], table_C[gene_row,5]),
#                        beside = true,
#                        ylim = c(0,max(table_C[gene_row,1],table_C[gene_row,2],table_C[gene_row,3], table_C[gene_row,4], table_C[gene_row,5])),
#                        cex.names = .75, names=c('0H', '2H','2H+JQ1','24H', '24H+JQ1') ,
#                        main = c(gene,cluster, "Cumulative"),
#                        ylab = "Average Signal",
#                        border = "black", axes = TRUE)
#     
# 
#   }
#   dev.off()
#   
#       
#   
# }


