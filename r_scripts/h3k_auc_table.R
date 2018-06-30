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

log_mat=cbind(c_matrix,log2(c_matrix[,2]/c_matrix[,1]),log2(c_matrix[,4]/c_matrix[,1]),log2(c_matrix[,4]/c_matrix[,2]))

colnames(log_mat)=c('H3K_0H','H3K_2H+PDGF','H3K_2H+PDGF+JQ1','H3K_24H+PDGF','H3K_2H+PDGF+JQ1','Log2(H3K_2H/H3K_0H)','Log2(H3K_24H/H3K_0H)','Log2(H3K_24H/H3K_2H)')

log_mat=round(log_mat,digits=2)
write.table(log_mat, file = "/storage/cylin/grail/projects/rasmc_all/tables/h3k27ac_gene_level_cumulative_signal.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")