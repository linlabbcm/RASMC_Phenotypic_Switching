a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4/RASMC_BRD4_UNSTIM_REP1_0_STITCH_-_JQ1/RASMC_BRD4_UNSTIM_REP1_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4/RASMC_BRD4_PDGF_2H_REP2_0_STITCH_-_JQ1/RASMC_BRD4_PDGF_2H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4/RASMC_BRD4_PDGF_24H_REP2_0_STITCH_-_JQ1/RASMC_BRD4_PDGF_24H_REP2_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)


a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE


promoterTable = cbind(a0H_GENE_TABLE[,2],a2H_GENE_TABLE[,2],a24H_GENE_TABLE[,2])
#print(promoterTable[1:10,])
enhancerTable = cbind(a0H_GENE_TABLE[,3],a2H_GENE_TABLE[,3],a24H_GENE_TABLE[,3])
#print(enhancerTable[1:10,])
cumulativeTable = cbind((promoterTable[,1]+enhancerTable[,1]),(promoterTable[,2]+enhancerTable[,2]),(promoterTable[,3]+enhancerTable[,3]))
#print(cumulativeTable[1:10,])

p_matrix = as.matrix(promoterTable)
e_matrix = as.matrix(enhancerTable)
c_matrix = as.matrix(cumulativeTable)


rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE

clusterCgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_4_genes_list.txt",header = FALSE)

clusterAgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_5_genes_list.txt",header = FALSE)

clusterBgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_1_genes_list.txt",header = FALSE)

clusterDgenes = read.delim("/storage/cylin/grail/projects/rasmc_all/tables/cluster_6_genes_list.txt",header = FALSE)


clusterArows = c()
clusterBrows = c()
clusterCrows = c()
clusterDrows = c()

 
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
for(i in 1:147){
  clusterDrow = which(a0h_matrix[,1] == clusterDgenes[i,1])
  clusterDrows = append(clusterDrows,clusterDrow)
}



A=clusterArows
B=clusterBrows
C=clusterCrows
D=clusterDrows
cluster_name1='Cluster A'
cluster_name2='Cluster B'
cluster_name3='Cluster C'
cluster_name4='Cluster D'



pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/BRD4_CLUSTER_A_B_C_D_promoter_enhancer_boxplots_old_reps.pdf')
par(mfrow=c(2,2))

############################################################################################################
table_PA <- p_matrix[clusterArows,]
boxplot(table_PA[,1], table_PA[,2], table_PA[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H', '24H'), ylab= c("Average Signal"), main = c('Promoter','Cluster A'))

table_PB <- p_matrix[clusterBrows,]
boxplot(table_PB[,1], table_PB[,2], table_PB[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H', '24H'),  ylab= c("Average Signal"), main = c('Promoter','Cluster B'))

table_PC <- p_matrix[clusterCrows,]
boxplot(table_PC[,1], table_PC[,2], table_PC[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H', '24H'), ylab= c("Average Signal"), main = c('Promoter','Cluster C'))

table_PD <- p_matrix[clusterDrows,]
boxplot(table_PD[,1], table_PD[,2], table_PD[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H','24H'), ylab= c("Average Signal"), main = c('Promoter','Cluster D'))

#########################################################################################################
table_EA <- e_matrix[clusterArows,]
boxplot(table_EA[,1], table_EA[,2], table_EA[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster A'))

table_EB <- e_matrix[clusterBrows,]
boxplot(table_EB[,1], table_EB[,2], table_EB[,3], cex=0, ylim=c(0,20000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster B'))

table_EC <- e_matrix[clusterCrows,]
boxplot(table_EC[,1], table_EC[,2], table_EC[,3], cex=0, ylim=c(0,20000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster C'))
table_ED <- e_matrix[clusterDrows,]
boxplot(table_ED[,1], table_ED[,2], table_ED[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H','24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster D'))


  ###########################################################################################################3
table_CA <- c_matrix[clusterArows,]
boxplot(table_CA[,1], table_CA[,2], table_CA[,3],cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster A'))

table_CB <- c_matrix[clusterBrows,]
boxplot(table_CB[,1], table_CB[,2], table_CB[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster B'))
table_CC <- c_matrix[clusterCrows,]
boxplot(table_CC[,1], table_CC[,2], table_CC[,3],cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster C'))
table_CD <- c_matrix[clusterDrows,]
boxplot(table_CD[,1], table_CD[,2], table_CD[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster D'))



dev.off()







pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/BRD4/BRD4_clusters_A_B_C_D_0h_normalized_enhancer_promoter_cumulative_boxplots_old_reps.pdf')
par(mfrow=c(2,2))



############################################################################################################################################

table_PA <- p_matrix[clusterArows,]
boxplot(log2(table_PA[,2]/table_PA[,1]), log2(table_PA[,3]/table_PA[,1]), cex=0, ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster A'))
table_PB <- p_matrix[clusterBrows,]
boxplot(log2(table_PB[,2]/table_PB[,1]), log2(table_PB[,3]/table_PB[,1]), cex=0,ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster B'))
table_PC <- p_matrix[clusterCrows,]
boxplot(log2(table_PC[,2]/table_PC[,1]), log2(table_PC[,3]/table_PC[,1]), cex=0, ylim=c(-15,10),names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster C'))
table_PD <- p_matrix[clusterDrows,]
boxplot(log2(table_PD[,2]/table_PD[,1]), log2(table_PD[,3]/table_PD[,1]), cex=0, ylim=c(-10,10), names = c( '2H','24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster D'))


  #############################################################################################################################################
table_EA <- e_matrix[clusterArows,]
boxplot(log2(table_EA[,2]/table_EA[,1]), log2(table_EA[,3]/table_EA[,1]), cex=0, ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster A'))
table_EB <- e_matrix[clusterBrows,]
boxplot(log2(table_EB[,2]/table_EB[,1]), log2(table_EB[,3]/table_EB[,1]), cex=0, ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster B'))
table_EC <- e_matrix[clusterCrows,]
boxplot(log2(table_EC[,2]/table_EC[,1]), log2(table_EC[,3]/table_EC[,1]), cex=0, ylim=c(-15,15), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster C'))
table_ED <- e_matrix[clusterDrows,]
boxplot(log2(table_ED[,2]/table_ED[,1]), log2(table_ED[,3]/table_ED[,1]), cex=0, ylim=c(-15,15), names = c( '2H','24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster D'))

  ############################################################################################################################################

table_CA <- c_matrix[clusterArows,]
boxplot(log2(table_CA[,2]/table_CA[,1]), log2(table_CA[,3]/table_CA[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster A'))
table_CB <- c_matrix[clusterBrows,]
boxplot(log2(table_CB[,2]/table_CB[,1]), log2(table_CB[,3]/table_CB[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster B'))
table_CC <- c_matrix[clusterCrows,]
boxplot(log2(table_CC[,2]/table_CC[,1]), log2(table_CC[,3]/table_CC[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster C'))
table_CD <- c_matrix[clusterDrows,]
boxplot(log2(table_CD[,2]/table_CD[,1]), log2(table_CD[,3]/table_CD[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster D'))

dev.off()


###################################################################
#############################BOXPOT T-TEST#########################
###################################################################
# t.test(table_CA[,1],table_CA[,2])
# t.test(table_CA[,1],table_CA[,3])
# t.test(table_CA[,2],table_CA[,3])
# 
# t.test(table_CB[,1],table_CB[,2])
# t.test(table_CB[,1],table_CB[,3])
# t.test(table_CB[,2],table_CB[,3])
# 
# t.test(table_CC[,1],table_CC[,2])
# t.test(table_CC[,1],table_CC[,3])
# t.test(table_CC[,2],table_CC[,3])



########################################################
#########################PIE CHARTS#####################
#########################################################

clusterAproActivity_0 = sum(table_PA[,1])
clusterAproActivity_2 = sum(table_PA[,2])
clusterAproActivity_24 = sum(table_PA[,3])

clusterAenhActivity_0 = sum(table_EA[,1])
clusterAenhActivity_2 = sum(table_EA[,2])
clusterAenhActivity_24 = sum(table_EA[,3])

A0p = clusterAproActivity_0/(sum(table_CA[,1]))
A0e = clusterAenhActivity_0/(sum(table_CA[,1]))

A2p = clusterAproActivity_2/(sum(table_CA[,2]))
A2e = clusterAenhActivity_2/(sum(table_CA[,2]))

A24p = clusterAproActivity_24/(sum(table_CA[,3]))
A24e = clusterAenhActivity_24/(sum(table_CA[,3]))

# ###########################################################
clusterBproActivity_0 = sum(table_PB[,1])
clusterBproActivity_2 = sum(table_PB[,2])
clusterBproActivity_24 = sum(table_PB[,3])

clusterBenhActivity_0 = sum(table_EB[,1])
clusterBenhActivity_2 = sum(table_EB[,2])
clusterBenhActivity_24 = sum(table_EB[,3])

B0p = clusterBproActivity_0/(sum(table_CB[,1]))
B0e = clusterBenhActivity_0/(sum(table_CB[,1]))

B2p = clusterBproActivity_2/(sum(table_CB[,2]))
B2e = clusterBenhActivity_2/(sum(table_CB[,2]))

B24p = clusterBproActivity_24/(sum(table_CB[,3]))
B24e = clusterBenhActivity_24/(sum(table_CB[,3]))

# ############################################################
clusterCproActivity_0 = sum(table_PC[,1])
clusterCproActivity_2 = sum(table_PC[,2])
clusterCproActivity_24 = sum(table_PC[,3])

clusterCenhActivity_0 = sum(table_EC[,1])
clusterCenhActivity_2 = sum(table_EC[,2])
clusterCenhActivity_24 = sum(table_EC[,3])

C0p = clusterCproActivity_0/(sum(table_CC[,1]))
C0e = clusterCenhActivity_0/(sum(table_CC[,1]))

C2p = clusterCproActivity_2/(sum(table_CC[,2]))
C2e = clusterCenhActivity_2/(sum(table_CC[,2]))

C24p = clusterCproActivity_24/(sum(table_CC[,3]))
C24e = clusterCenhActivity_24/(sum(table_CC[,3]))

# ##############################################################

clusterDproActivity_0 = sum(table_PD[,1])
clusterDproActivity_2 = sum(table_PD[,2])
clusterDproActivity_24 = sum(table_PD[,3])

clusterDenhActivity_0 = sum(table_ED[,1])
clusterDenhActivity_2 = sum(table_ED[,2])
clusterDenhActivity_24 = sum(table_ED[,3])

D0p = clusterDproActivity_0/(sum(table_CD[,1]))
D0e = clusterDenhActivity_0/(sum(table_CD[,1]))

D2p = clusterDproActivity_2/(sum(table_CD[,2]))
D2e = clusterDenhActivity_2/(sum(table_CD[,2]))

D24p = clusterDproActivity_24/(sum(table_CD[,3]))
D24e = clusterDenhActivity_24/(sum(table_CD[,3]))

# ##############################################################
table_P = p_matrix
table_E = e_matrix
table_C = c_matrix

proActivity_0 = sum(table_P[,1])
proActivity_2 = sum(table_P[,2])
proActivity_24 = sum(table_P[,3])

enhActivity_0 = sum(table_E[,1])
enhActivity_2 = sum(table_E[,2])
enhActivity_24 = sum(table_E[,3])

All0p = proActivity_0/(sum(table_C[,1]))
All0e = enhActivity_0/(sum(table_C[,1]))

print("All0p")
print(All0p)
print("All0e")
print(All0e)

All2p = proActivity_2/(sum(table_C[,2]))
All2e = enhActivity_2/(sum(table_C[,2]))

print("All2p")
print(All2p)
print("All2e")
print(All2e)


All24p = proActivity_24/(sum(table_C[,3]))
All24e = enhActivity_24/(sum(table_C[,3]))

print("All24p")
print(All24p)
print("All24e")
print(All24e)

###########################################################

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/ALL_ACTIVE_GENES_BRD4_ENH_vs_PRO_PIE_CHARTS_old_reps.pdf')
par(mfrow=c(3,1))


slices <- c(All0p,All0e)
lbls <- c(All0p,All0e)
colors = c('grey','black')
pie(slices,labels=lbls,col = colors,
    main=c("All Active Genes","Promoter(grey) vs Enhancer(black) Activity","0H BRD4"))

slices <- c(All2p,All2e)
lbls <- c(All2p,All2e)
colors = c('grey','black')
pie(slices,labels=lbls,col = colors,
    main=c("All Active Genes","Promoter(grey) vs Enhancer(black) Activity","2H BRD4"))

slices <- c(All24p,All24e)
lbls <- c(All24p,All24e)
colors = c('grey','black')
pie(slices,labels=lbls,col = colors,
    main=c("All Active Genes","Promoter(grey) vs Enhancer(black) Activity","24H BRD4"))

dev.off()

###########################################################
#
pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/CLUSTERS_BRD4_ENH_vs_PRO_PIE_CHARTS_old_reps.pdf')
par(mfrow=c(3,1))


slices <- c(A0p,A0e)
lbls <- c(A0p,A0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster A","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(A2p,A2e)
lbls <- c(A2p,A2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster A","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(A24p,A24e)
lbls <- c(A24p,A24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster A","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
###########################################################################
slices <- c(B0p,B0e)
lbls <- c(B0p,B0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster B","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(B2p,B2e)
lbls <- c(B2p,B2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster B","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(B24p,B24e)
lbls <- c(B24p,B24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster B","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
###########################################################################
slices <- c(C0p,C0e)
lbls <- c(C0p,C0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
      main=c("Cluster C","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(C2p,C2e)
lbls <- c(C2p,C2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster C","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(C24p,C24e)
lbls <- c(C24p,C24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster C","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
#############################################################################
slices <- c(D0p,D0e)
lbls <- c(D0p,D0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster D","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(D2p,D2e)
lbls <- c(D2p,D2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster D","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(D24p,D24e)
lbls <- c(D24p,D24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster D","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
###########################################################################


dev.off()




##########################################################################
#####################BRD4 BARGRAPHS#######################################
##########################################################################
# 
# gene_name = c('E2f4','Thra','Rara','Nr1d1','Osr1')
# 
# pdf('/storage/cylin/grail/projects/rasmc_all/figures/Figure_2_Brd4_BARGRAPHS.pdf')
# par(mfrow=c(3,2))
# 
# for(gene in gene_name){ 
#   gene_row = which(rownames(c_matrix) == gene)
#   
#   print(gene_row)  
#   
#   if(length(gene_row) != 0){ 
#     plot_0H = c()
#     plot_2H =  c()
#     plot_24H = c()
#     
#     
#     plot_0H = c((a0H_GENE_TABLE[gene_row,2]+a0H_GENE_TABLE[gene_row,3]),(b0H_GENE_TABLE[gene_row,2]+b0H_GENE_TABLE[gene_row,3]))
#     plot_2H =  c((a2H_GENE_TABLE[gene_row,2]+a2H_GENE_TABLE[gene_row,3]),(b2H_GENE_TABLE[gene_row,2]+b2H_GENE_TABLE[gene_row,3]))
#     plot_24H = c((a24H_GENE_TABLE[gene_row,2]+a24H_GENE_TABLE[gene_row,3]),(b24H_GENE_TABLE[gene_row,2]+b24H_GENE_TABLE[gene_row,3]))
#     
#     mean_0H = mean(plot_0H)
#     mean_2H = mean(plot_2H)
#     mean_24H = mean(plot_24H)
#     
#     sd_0H = sd(plot_0H)
#     sd_2H = sd(plot_2H)
#     sd_24H = sd(plot_24H)
#     
#     se_0H = sd_0H/(sqrt(length(plot_0H)))
#     se_2H = sd_2H/(sqrt(length(plot_2H)))
#     se_24H = sd_24H/(sqrt(length(plot_24H)))
#     
#     maxy=max(mean_0H,mean_2H, mean_24H) + 5*sd_24H 
#     print(maxy) 
#     
#     centers <- barplot(height = c(mean_0H,mean_2H, mean_24H),
#                        beside = true, las = 2,
#                        ylim = c(0, maxy),
#                        cex.names = 0.75, 
#                        main = gene,
#                        ylab = "Brd4 Signal",
#                        border = "black", axes = TRUE)
#     
#     
#     # 45 degree string rotation
#     text(x = centers, y = par("usr")[3] - 1, srt = 45,
#          adj = 1, labels = c('0H', '2H','24H'), xpd = TRUE)
#     
#     segments(centers, c(mean_0H,mean_2H, mean_24H) - c(se_0H,se_2H, se_24H) * 2, centers,
#              c(mean_0H,mean_2H, mean_24H) + c(se_0H,se_2H,se_24H) * 2, lwd = 1.5)
#     
#     arrows(centers, c(mean_0H,mean_2H, mean_24H) - c(se_0H,se_2H,se_24H) * 2, centers,
#            c(mean_0H,mean_2H, mean_24H) + c(se_0H,se_2H,se_24H) * 2, lwd = 1.5, angle = 90,
#            code = 3, length = 0.05)
#   }
# }
# 
# dev.off()


b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4/RASMC_BRD4_UNSTIM_NEW_0_STITCH_-_JQ1/RASMC_BRD4_UNSTIM_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4/RASMC_BRD4_PDGF_24H_NEW_0_STITCH_-_JQ1/RASMC_BRD4_PDGF_24H_NEW_0_STITCH_-_JQ1_GENE_TABLE.txt', sep = '\t', header = TRUE)

promoterTable = cbind(b0H_GENE_TABLE[,2],b2H_GENE_TABLE[,2],b24H_GENE_TABLE[,2])
#print(promoterTable[1:10,])
enhancerTable = cbind(b0H_GENE_TABLE[,3],b2H_GENE_TABLE[,3],b24H_GENE_TABLE[,3])
#print(enhancerTable[1:10,])
cumulativeTable = cbind((promoterTable[,1]+enhancerTable[,1]),(promoterTable[,2]+enhancerTable[,2]),(promoterTable[,3]+enhancerTable[,3]))
#print(cumulativeTable[1:10,])

p_matrix = as.matrix(promoterTable)
e_matrix = as.matrix(enhancerTable)
c_matrix = as.matrix(cumulativeTable)


rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE

clusterArows = c()
clusterBrows = c()
clusterCrows = c()
clusterDrows = c()
# 
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
for(i in 1:147){
  clusterDrow = which(a0h_matrix[,1] == clusterDgenes[i,1])
  clusterDrows = append(clusterDrows,clusterDrow)
}



A=clusterArows
B=clusterBrows
C=clusterCrows
D=clusterDrows
cluster_name1='Cluster A'
cluster_name2='Cluster B'
cluster_name3='Cluster C'
cluster_name4='Cluster D'



pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/BRD4_CLUSTER_A_B_C_D_promoter_enhancer_boxplots_new_reps.pdf')
par(mfrow=c(2,2))

############################################################################################################
table_PA <- p_matrix[clusterArows,]
boxplot(table_PA[,1], table_PA[,2], table_PA[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H', '24H'), ylab= c("Average Signal"), main = c('Promoter','Cluster A'))

table_PB <- p_matrix[clusterBrows,]
boxplot(table_PB[,1], table_PB[,2], table_PB[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H', '24H'),  ylab= c("Average Signal"), main = c('Promoter','Cluster B'))

table_PC <- p_matrix[clusterCrows,]
boxplot(table_PC[,1], table_PC[,2], table_PC[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H', '24H'), ylab= c("Average Signal"), main = c('Promoter','Cluster C'))

table_PD <- p_matrix[clusterDrows,]
boxplot(table_PD[,1], table_PD[,2], table_PD[,3], cex=0, ylim=c(0,3000), names = c('0H', '2H','24H'), ylab= c("Average Signal"), main = c('Promoter','Cluster D'))

#########################################################################################################
table_EA <- e_matrix[clusterArows,]
boxplot(table_EA[,1], table_EA[,2], table_EA[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster A'))

table_EB <- e_matrix[clusterBrows,]
boxplot(table_EB[,1], table_EB[,2], table_EB[,3], cex=0, ylim=c(0,20000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster B'))

table_EC <- e_matrix[clusterCrows,]
boxplot(table_EC[,1], table_EC[,2], table_EC[,3], cex=0, ylim=c(0,20000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster C'))
table_ED <- e_matrix[clusterDrows,]
boxplot(table_ED[,1], table_ED[,2], table_ED[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H','24H'),
          ylab= c("Average Signal"), main = c('Enhancer','Cluster D'))


  ###########################################################################################################3
table_CA <- c_matrix[clusterArows,]
boxplot(table_CA[,1], table_CA[,2], table_CA[,3],cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster A'))

table_CB <- c_matrix[clusterBrows,]
boxplot(table_CB[,1], table_CB[,2], table_CB[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster B'))
table_CC <- c_matrix[clusterCrows,]
boxplot(table_CC[,1], table_CC[,2], table_CC[,3],cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster C'))
table_CD <- c_matrix[clusterDrows,]
boxplot(table_CD[,1], table_CD[,2], table_CD[,3], cex=0, ylim=c(0,25000), names = c('0H', '2H', '24H'),
          ylab= c("Average Signal"), main = c('Cumulative','Cluster D'))



dev.off()







pdf(file = '/storage/cylin/grail/projects/rasmc_all/figures/BRD4/BRD4_clusters_A_B_C_D_0h_normalized_enhancer_promoter_cumulative_boxplots_new_reps.pdf')
par(mfrow=c(2,2))


############################################################################################################################################

table_PA <- p_matrix[clusterArows,]
boxplot(log2(table_PA[,2]/table_PA[,1]), log2(table_PA[,3]/table_PA[,1]), cex=0, ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster A'))
table_PB <- p_matrix[clusterBrows,]
boxplot(log2(table_PB[,2]/table_PB[,1]), log2(table_PB[,3]/table_PB[,1]), cex=0,ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster B'))
table_PC <- p_matrix[clusterCrows,]
boxplot(log2(table_PC[,2]/table_PC[,1]), log2(table_PC[,3]/table_PC[,1]), cex=0, ylim=c(-15,10),names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster C'))
table_PD <- p_matrix[clusterDrows,]
boxplot(log2(table_PD[,2]/table_PD[,1]), log2(table_PD[,3]/table_PD[,1]), cex=0, ylim=c(-10,10), names = c( '2H','24H'),
          ylab= c("0H Normalized Signal"), main = c('Promoter','Cluster D'))


  #############################################################################################################################################
table_EA <- e_matrix[clusterArows,]
boxplot(log2(table_EA[,2]/table_EA[,1]), log2(table_EA[,3]/table_EA[,1]), cex=0, ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster A'))
table_EB <- e_matrix[clusterBrows,]
boxplot(log2(table_EB[,2]/table_EB[,1]), log2(table_EB[,3]/table_EB[,1]), cex=0, ylim=c(-15,10), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster B'))
table_EC <- e_matrix[clusterCrows,]
boxplot(log2(table_EC[,2]/table_EC[,1]), log2(table_EC[,3]/table_EC[,1]), cex=0, ylim=c(-15,15), names = c( '2H', '24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster C'))
table_ED <- e_matrix[clusterDrows,]
boxplot(log2(table_ED[,2]/table_ED[,1]), log2(table_ED[,3]/table_ED[,1]), cex=0, ylim=c(-15,15), names = c( '2H','24H'),
          ylab= c("0H Normalized Signal"), main = c('Enhancer','Cluster D'))

  ############################################################################################################################################

table_CA <- c_matrix[clusterArows,]
boxplot(log2(table_CA[,2]/table_CA[,1]), log2(table_CA[,3]/table_CA[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster A'))
table_CB <- c_matrix[clusterBrows,]
boxplot(log2(table_CB[,2]/table_CB[,1]), log2(table_CB[,3]/table_CB[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster B'))
table_CC <- c_matrix[clusterCrows,]
boxplot(log2(table_CC[,2]/table_CC[,1]), log2(table_CC[,3]/table_CC[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster C'))
table_CD <- c_matrix[clusterDrows,]
boxplot(log2(table_CD[,2]/table_CD[,1]), log2(table_CD[,3]/table_CD[,1]), cex=0, ylim=c(-10,10), names = c( '2H', '24H'),
          ylab= c("Log2 Fold Change 0H Normalized Signal"), main = c('Cumulative','Cluster D'))

dev.off()


###################################################################
#############################BOXPOT T-TEST#########################
###################################################################
# t.test(table_CA[,1],table_CA[,2])
# t.test(table_CA[,1],table_CA[,3])
# t.test(table_CA[,2],table_CA[,3])
# 
# t.test(table_CB[,1],table_CB[,2])
# t.test(table_CB[,1],table_CB[,3])
# t.test(table_CB[,2],table_CB[,3])
# 
# t.test(table_CC[,1],table_CC[,2])
# t.test(table_CC[,1],table_CC[,3])
# t.test(table_CC[,2],table_CC[,3])


########################################################
#########################PIE CHARTS#####################
#########################################################

clusterAproActivity_0 = sum(table_PA[,1])
clusterAproActivity_2 = sum(table_PA[,2])
clusterAproActivity_24 = sum(table_PA[,3])

clusterAenhActivity_0 = sum(table_EA[,1])
clusterAenhActivity_2 = sum(table_EA[,2])
clusterAenhActivity_24 = sum(table_EA[,3])

A0p = clusterAproActivity_0/(sum(table_CA[,1]))
A0e = clusterAenhActivity_0/(sum(table_CA[,1]))

A2p = clusterAproActivity_2/(sum(table_CA[,2]))
A2e = clusterAenhActivity_2/(sum(table_CA[,2]))

A24p = clusterAproActivity_24/(sum(table_CA[,3]))
A24e = clusterAenhActivity_24/(sum(table_CA[,3]))

# ###########################################################
clusterBproActivity_0 = sum(table_PB[,1])
clusterBproActivity_2 = sum(table_PB[,2])
clusterBproActivity_24 = sum(table_PB[,3])

clusterBenhActivity_0 = sum(table_EB[,1])
clusterBenhActivity_2 = sum(table_EB[,2])
clusterBenhActivity_24 = sum(table_EB[,3])

B0p = clusterBproActivity_0/(sum(table_CB[,1]))
B0e = clusterBenhActivity_0/(sum(table_CB[,1]))

B2p = clusterBproActivity_2/(sum(table_CB[,2]))
B2e = clusterBenhActivity_2/(sum(table_CB[,2]))

B24p = clusterBproActivity_24/(sum(table_CB[,3]))
B24e = clusterBenhActivity_24/(sum(table_CB[,3]))

# ############################################################
clusterCproActivity_0 = sum(table_PC[,1])
clusterCproActivity_2 = sum(table_PC[,2])
clusterCproActivity_24 = sum(table_PC[,3])

clusterCenhActivity_0 = sum(table_EC[,1])
clusterCenhActivity_2 = sum(table_EC[,2])
clusterCenhActivity_24 = sum(table_EC[,3])

C0p = clusterCproActivity_0/(sum(table_CC[,1]))
C0e = clusterCenhActivity_0/(sum(table_CC[,1]))

C2p = clusterCproActivity_2/(sum(table_CC[,2]))
C2e = clusterCenhActivity_2/(sum(table_CC[,2]))

C24p = clusterCproActivity_24/(sum(table_CC[,3]))
C24e = clusterCenhActivity_24/(sum(table_CC[,3]))

# ##############################################################

clusterDproActivity_0 = sum(table_PD[,1])
clusterDproActivity_2 = sum(table_PD[,2])
clusterDproActivity_24 = sum(table_PD[,3])

clusterDenhActivity_0 = sum(table_ED[,1])
clusterDenhActivity_2 = sum(table_ED[,2])
clusterDenhActivity_24 = sum(table_ED[,3])

D0p = clusterDproActivity_0/(sum(table_CD[,1]))
D0e = clusterDenhActivity_0/(sum(table_CD[,1]))

D2p = clusterDproActivity_2/(sum(table_CD[,2]))
D2e = clusterDenhActivity_2/(sum(table_CD[,2]))

D24p = clusterDproActivity_24/(sum(table_CD[,3]))
D24e = clusterDenhActivity_24/(sum(table_CD[,3]))

# ##############################################################
table_P = p_matrix
table_E = e_matrix
table_C = c_matrix

proActivity_0 = sum(table_P[,1])
proActivity_2 = sum(table_P[,2])
proActivity_24 = sum(table_P[,3])

enhActivity_0 = sum(table_E[,1])
enhActivity_2 = sum(table_E[,2])
enhActivity_24 = sum(table_E[,3])

All0p = proActivity_0/(sum(table_C[,1]))
All0e = enhActivity_0/(sum(table_C[,1]))

print("All0p")
print(All0p)
print("All0e")
print(All0e)

All2p = proActivity_2/(sum(table_C[,2]))
All2e = enhActivity_2/(sum(table_C[,2]))

print("All2p")
print(All2p)
print("All2e")
print(All2e)


All24p = proActivity_24/(sum(table_C[,3]))
All24e = enhActivity_24/(sum(table_C[,3]))

print("All24p")
print(All24p)
print("All24e")
print(All24e)

###########################################################

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/ALL_ACTIVE_GENES_BRD4_ENH_vs_PRO_PIE_CHARTS_new_reps.pdf')
par(mfrow=c(3,1))


slices <- c(All0p,All0e)
lbls <- c(All0p,All0e)
colors = c('grey','black')
pie(slices,labels=lbls,col = colors,
    main=c("All Active Genes","Promoter(grey) vs Enhancer(black) Activity","0H BRD4"))

slices <- c(All2p,All2e)
lbls <- c(All2p,All2e)
colors = c('grey','black')
pie(slices,labels=lbls,col = colors,
    main=c("All Active Genes","Promoter(grey) vs Enhancer(black) Activity","2H BRD4"))

slices <- c(All24p,All24e)
lbls <- c(All24p,All24e)
colors = c('grey','black')
pie(slices,labels=lbls,col = colors,
    main=c("All Active Genes","Promoter(grey) vs Enhancer(black) Activity","24H BRD4"))

dev.off()

###########################################################
#
pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/CLUSTERS_BRD4_ENH_vs_PRO_PIE_CHARTS_new_reps.pdf')
par(mfrow=c(3,1))


slices <- c(A0p,A0e)
lbls <- c(A0p,A0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster A","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(A2p,A2e)
lbls <- c(A2p,A2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster A","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(A24p,A24e)
lbls <- c(A24p,A24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster A","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
###########################################################################
slices <- c(B0p,B0e)
lbls <- c(B0p,B0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster B","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(B2p,B2e)
lbls <- c(B2p,B2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster B","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(B24p,B24e)
lbls <- c(B24p,B24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster B","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
###########################################################################
slices <- c(C0p,C0e)
lbls <- c(C0p,C0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
      main=c("Cluster C","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(C2p,C2e)
lbls <- c(C2p,C2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster C","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(C24p,C24e)
lbls <- c(C24p,C24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster C","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
#############################################################################
slices <- c(D0p,D0e)
lbls <- c(D0p,D0e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster D","Promoter(red) vs Enhancer(blue) Activity","0H BRD4"))

slices <- c(D2p,D2e)
lbls <- c(D2p,D2e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster D","Promoter(red) vs Enhancer(blue) Activity","2H BRD4"))

slices <- c(D24p,D24e)
lbls <- c(D24p,D24e)
colors = c('red','blue')
pie(slices,labels=lbls,col = colors,
    main=c("Cluster D","Promoter(red) vs Enhancer(blue) Activity","24H BRD4"))
###########################################################################


dev.off()

