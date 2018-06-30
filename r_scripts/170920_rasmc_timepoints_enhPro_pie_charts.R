
a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_UNSTIM_h3k/RASMC_BRD4_UNSTIM_REP1_h3k_unstim/RASMC_BRD4_UNSTIM_REP1_h3k_unstim_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_UNSTIM_h3k/RASMC_BRD4_PDGF_2H_REP2_0_h3k_unstim/RASMC_BRD4_PDGF_2H_REP2_0_h3k_unstim_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_UNSTIM_h3k/RASMC_BRD4_PDGF_24H_REP2_h3k_unstim/RASMC_BRD4_PDGF_24H_REP2_h3k_unstim_GENE_TABLE.txt', sep = '\t', header = TRUE)
b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_UNSTIM_h3k/RASMC_BRD4_UNSTIM_NEW_h3k_unstim/RASMC_BRD4_UNSTIM_NEW_h3k_unstim_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_UNSTIM_h3k/RASMC_BRD4_PDGF_2H_NEW_h3k_unstim/RASMC_BRD4_PDGF_2H_NEW_h3k_unstim_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_UNSTIM_h3k/RASMC_BRD4_PDGF_24H_h3k_unstim/RASMC_BRD4_PDGF_24H_h3k_unstim_GENE_TABLE.txt', sep = '\t', header = TRUE)


a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE


promoterTable = cbind((a0H_GENE_TABLE[,2]+b0H_GENE_TABLE[,2])/2,(a2H_GENE_TABLE[,2]+b2H_GENE_TABLE[,2])/2,(a24H_GENE_TABLE[,2]+b24H_GENE_TABLE[,2])/2)
#print(promoterTable[1:10,])
enhancerTable = cbind((a0H_GENE_TABLE[,3]+b0H_GENE_TABLE[,3])/2,(a2H_GENE_TABLE[,3]+b2H_GENE_TABLE[,3])/2,(a24H_GENE_TABLE[,3]+b24H_GENE_TABLE[,3])/2)
#print(enhancerTable[1:10,])
cumulativeTable = cbind((promoterTable[,1]+enhancerTable[,1]),(promoterTable[,2]+enhancerTable[,2]),(promoterTable[,3]+enhancerTable[,3]))
#print(cumulativeTable[1:10,])

p_matrix = as.matrix(promoterTable)
e_matrix = as.matrix(enhancerTable)
c_matrix = as.matrix(cumulativeTable)


rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE






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
# 
###########################################################

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/ALL_ACTIVE_GENES_BRD4_ENH_vs_PRO_PIE_CHARTS_h3k_unstim.pdf')
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







a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_2H_h3k/RASMC_BRD4_UNSTIM_REP1_h3k_2h/RASMC_BRD4_UNSTIM_REP1_h3k_2h_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_2H_h3k/RASMC_BRD4_PDGF_2H_REP2_0_h3k_2h/RASMC_BRD4_PDGF_2H_REP2_0_h3k_2h_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_2H_h3k/RASMC_BRD4_PDGF_24H_REP2_h3k_2h/RASMC_BRD4_PDGF_24H_REP2_h3k_2h_GENE_TABLE.txt', sep = '\t', header = TRUE)
b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_2H_h3k/RASMC_BRD4_UNSTIM_NEW_h3k_2h/RASMC_BRD4_UNSTIM_NEW_h3k_2h_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_2H_h3k/RASMC_BRD4_PDGF_2H_NEW_h3k_2h/RASMC_BRD4_PDGF_2H_NEW_h3k_2h_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_2H_h3k/RASMC_BRD4_PDGF_24H_h3k_2h/RASMC_BRD4_PDGF_24H_h3k_2h_GENE_TABLE.txt', sep = '\t', header = TRUE)


a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE


promoterTable = cbind((a0H_GENE_TABLE[,2]+b0H_GENE_TABLE[,2])/2,(a2H_GENE_TABLE[,2]+b2H_GENE_TABLE[,2])/2,(a24H_GENE_TABLE[,2]+b24H_GENE_TABLE[,2])/2)
#print(promoterTable[1:10,])
enhancerTable = cbind((a0H_GENE_TABLE[,3]+b0H_GENE_TABLE[,3])/2,(a2H_GENE_TABLE[,3]+b2H_GENE_TABLE[,3])/2,(a24H_GENE_TABLE[,3]+b24H_GENE_TABLE[,3])/2)
#print(enhancerTable[1:10,])
cumulativeTable = cbind((promoterTable[,1]+enhancerTable[,1]),(promoterTable[,2]+enhancerTable[,2]),(promoterTable[,3]+enhancerTable[,3]))
#print(cumulativeTable[1:10,])

p_matrix = as.matrix(promoterTable)
e_matrix = as.matrix(enhancerTable)
c_matrix = as.matrix(cumulativeTable)


rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE






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
# 
###########################################################

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/ALL_ACTIVE_GENES_BRD4_ENH_vs_PRO_PIE_CHARTS_h3k_2h.pdf')
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










a0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_24H_h3k/RASMC_BRD4_UNSTIM_REP1_h3k_24h/RASMC_BRD4_UNSTIM_REP1_h3k_24h_GENE_TABLE.txt', sep = '\t', header = TRUE)
a2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_24H_h3k/RASMC_BRD4_PDGF_2H_REP2_0_h3k_24h/RASMC_BRD4_PDGF_2H_REP2_0_h3k_24h_GENE_TABLE.txt', sep = '\t', header = TRUE)
a24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_24H_h3k/RASMC_BRD4_PDGF_24H_REP2_h3k_24h/RASMC_BRD4_PDGF_24H_REP2_h3k_24h_GENE_TABLE.txt', sep = '\t', header = TRUE)
b0H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_24H_h3k/RASMC_BRD4_UNSTIM_NEW_h3k_24h/RASMC_BRD4_UNSTIM_NEW_h3k_24h_GENE_TABLE.txt', sep = '\t', header = TRUE)
b2H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_24H_h3k/RASMC_BRD4_PDGF_2H_NEW_h3k_24h/RASMC_BRD4_PDGF_2H_NEW_h3k_24h_GENE_TABLE.txt', sep = '\t', header = TRUE)
b24H_GENE_TABLE = read.delim('/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_24H_h3k/RASMC_BRD4_PDGF_24H_h3k_24h/RASMC_BRD4_PDGF_24H_h3k_24h_GENE_TABLE.txt', sep = '\t', header = TRUE)


a0h_matrix = as.matrix(a0H_GENE_TABLE)
genes = a0H_GENE_TABLE$GENE


promoterTable = cbind((a0H_GENE_TABLE[,2]+b0H_GENE_TABLE[,2])/2,(a2H_GENE_TABLE[,2]+b2H_GENE_TABLE[,2])/2,(a24H_GENE_TABLE[,2]+b24H_GENE_TABLE[,2])/2)
#print(promoterTable[1:10,])
enhancerTable = cbind((a0H_GENE_TABLE[,3]+b0H_GENE_TABLE[,3])/2,(a2H_GENE_TABLE[,3]+b2H_GENE_TABLE[,3])/2,(a24H_GENE_TABLE[,3]+b24H_GENE_TABLE[,3])/2)
#print(enhancerTable[1:10,])
cumulativeTable = cbind((promoterTable[,1]+enhancerTable[,1]),(promoterTable[,2]+enhancerTable[,2]),(promoterTable[,3]+enhancerTable[,3]))
#print(cumulativeTable[1:10,])

p_matrix = as.matrix(promoterTable)
e_matrix = as.matrix(enhancerTable)
c_matrix = as.matrix(cumulativeTable)


rownames(p_matrix) = a0H_GENE_TABLE$GENE
rownames(e_matrix) = a0H_GENE_TABLE$GENE
rownames(c_matrix) = a0H_GENE_TABLE$GENE






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
# 
###########################################################

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/BRD4/ALL_ACTIVE_GENES_BRD4_ENH_vs_PRO_PIE_CHARTS_h3k_24h.pdf')
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
