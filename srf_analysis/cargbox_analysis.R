# library(BSgenome.Rnorvegicus.UCSC.rn6)
# 
# 
# dnstr <- getSeq(genome, gr.peaks.rel)
# m.alpha.freq <- dinucleotideFrequency(dnstr)
# gc <- (m.alpha.freq[ ,"G"] + m.alpha.freq[ ,"C"])/(m.alpha.freq[ ,"G"] + m.alpha.freq[ ,"C"] +
#                                                      m.alpha.freq[ ,"A"] + m.alpha.freq[ ,"T"])
# gc <- round(gc, digit = 2) * 100
# df.gc <- data.frame(GC = gc)

freq_table=read.delim('C://Users/rhirsch/Documents/srf_analysis/carg_box_freq_table.txt')
bg_table=read.delim('C://Users/rhirsch/Documents/srf_analysis/all_motif_bgs.txt')
freq_mat = as.matrix(freq_table[,2:5])
rownames(freq_mat)=freq_table[,1]

bg_mat=as.matrix(bg_table[,2:17])
rownames(bg_mat)=bg_table[,1]


ex_val = (bg_mat[,6]*((bg_mat[,1]+bg_mat[,4]+bg_mat[,13]+bg_mat[,16])^3)*bg_mat[,11])*1000

log_vec = log2(freq_mat[,3]/(ex_val))

freq_mat=cbind(freq_mat,log_vec)

log_order = order(freq_mat[,5])

freq_mat[,4]=freq_mat[,4]*100

pdf(file='C://Users/rhirsch/Documents/srf_analysis/srf_enrichment_barplot.pdf')
barCenters <- barplot(freq_mat[log_order,5],las=2,main='Log2 FC actual vs expect Srf motif co-occupancy',ylim=c(0,8),names=rownames(freq_mat[log_order,]),ylab='Log2 FC actual vs expected')
segments(barCenters, freq_mat[log_order,5] - freq_mat[log_order,4], barCenters, freq_mat[log_order,5] + freq_mat[log_order,4], lwd = 1.5)
arrows(barCenters, freq_mat[log_order,5] - freq_mat[log_order,4], barCenters,freq_mat[log_order,5] + freq_mat[log_order,4],lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()
