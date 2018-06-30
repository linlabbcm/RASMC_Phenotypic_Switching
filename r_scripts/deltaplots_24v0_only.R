options(scipen=999)

#READ IN AND PLOT 2H vs 0H DELTA OUT TABLE

delta_table_2v0 = read.delim('C:/Users/rhirsch/Documents/am_ef/rasmc_h3k27ac_0_tss_am_ef_EDGE_DELTA_OUT_2V0.txt')

delta_matrix_2v0 = delta_table_2v0[,2:6]
rownames(delta_matrix_2v0) = delta_table_2v0[,1]
delta_matrix_2v0=as.matrix(delta_matrix_2v0)


#READ IN AND PLOT 24H vs 0H DELTA OUT TABLE

delta_table_24v0 = read.delim('C:/Users/rhirsch/Documents/am_ef/rasmc_h3k27ac_0_tss_am_ef_EDGE_DELTA_OUT_24V0.txt')

delta_matrix_24v0 = delta_table_24v0[,2:6]
rownames(delta_matrix_24v0) = delta_table_24v0[,1]
delta_matrix_24v0=as.matrix(delta_matrix_24v0)

#READ IN AND PLOT 24H vs 2H DELTA OUT TABLE

delta_table_24v2 = read.delim('C:/Users/rhirsch/Documents/am_ef/rasmc_h3k27ac_0_tss_am_ef_EDGE_DELTA_OUT_24V2.txt')

delta_matrix_24v2 = delta_table_24v2[,2:6]
rownames(delta_matrix_24v2) = delta_table_24v2[,1]
delta_matrix_24v2=as.matrix(delta_matrix_24v2)


#GRAB TF NAMES
tf_names=as.character(delta_table_24v0[,1])


out_order_24v0=order(delta_matrix_24v0[,3])

pdf(file='C:/Users/rhirsch/Documents/am_ef/out_degree_log2_fc_brd4_24v0_error_bars.pdf',width=18)

barCenters <- barplot(delta_matrix_24v0[out_order_24v0,3],las=2,main='Rank order by log2 FC Brd4 out degree 24H vs 0H',ylim=c(-.5,.5))
segments(barCenters, delta_matrix_24v0[out_order_24v0,3] - delta_matrix_24v0[out_order_24v0,5], barCenters,delta_matrix_24v0[out_order_24v0,3] + delta_matrix_24v0[out_order_24v0,5], lwd = 1.5)
arrows(barCenters, delta_matrix_24v0[out_order_24v0,3] - delta_matrix_24v0[out_order_24v0,5], barCenters,delta_matrix_24v0[out_order_24v0,3] + delta_matrix_24v0[out_order_24v0,5], lwd = 1.5, angle = 90,code = 3, length = 0.05)

dev.off()

out_order_2v0=order(delta_matrix_2v0[,3])

pdf(file='C:/Users/rhirsch/Documents/am_ef/out_degree_log2_fc_brd4_2v0_error_bars.pdf',width=18)

barCenters <- barplot(delta_matrix_2v0[out_order_2v0,3],las=2,main='Rank order by log2 FC Brd4 out degree 2H vs 0H',ylim=c(-.1,1.5))
segments(barCenters, delta_matrix_2v0[out_order_2v0,3] - delta_matrix_2v0[out_order_2v0,5], barCenters,delta_matrix_2v0[out_order_2v0,3] + delta_matrix_2v0[out_order_2v0,5], lwd = 1.5)
arrows(barCenters, delta_matrix_2v0[out_order_2v0,3] - delta_matrix_2v0[out_order_2v0,5], barCenters,delta_matrix_2v0[out_order_2v0,3] + delta_matrix_2v0[out_order_2v0,5], lwd = 1.5, angle = 90,code = 3, length = 0.05)

dev.off()

out_order_24v2=order(delta_matrix_24v2[,3])

pdf(file='C:/Users/rhirsch/Documents/am_ef/out_degree_log2_fc_brd4_24v2_error_bars.pdf',width=18)

barCenters <- barplot(delta_matrix_24v2[out_order_24v2,3],las=2,main='Rank order by log2 FC Brd4 out degree 24H vs 2H',ylim=c(-1.5,.1))
segments(barCenters, delta_matrix_24v2[out_order_24v2,3] - delta_matrix_24v2[out_order_24v2,5], barCenters,delta_matrix_24v2[out_order_24v2,3] + delta_matrix_24v2[out_order_24v2,5], lwd = 1.5)
arrows(barCenters, delta_matrix_24v2[out_order_24v2,3] - delta_matrix_24v2[out_order_24v2,5], barCenters,delta_matrix_24v2[out_order_24v2,3] + delta_matrix_24v2[out_order_24v2,5], lwd = 1.5, angle = 90,code = 3, length = 0.05)

dev.off()


########################################################
#################CROP FOR STRING TFS####################
########################################################

string_tfs=c('Dbp','Arntl','Bhlhe40','Id1','Tcf3','Tead3','Tead1','Pbx1','Meis1','Nfe2l2',
             'Maff','Mafg','Mafk','Thrb','Nr1d1','Rxrb','Thra','Esrra','Rarg','Nr2f6',
             'Rara','Rxra','Srebf1','Nr2f2','Egr2','Cebpd','Mnt','Myc','Yy1','Sp1','Smad6',
             'Smad3','Smad7','Cebpb','Rela','Nfkb2','Stat3','Foxo1','Irf9','Stat6',
             'Atf4','Jun','Ets1','Jund','Junb','Fosl1','Fosl2','Atf7','Irf3','Irf7')

string_tfs=sort(string_tfs)
rows=c()
for(tf in string_tfs){
  row = which(rownames(in_degree_Matrix)==toupper(tf))
  rows = c(rows,row)
}



rows=c()
for(tf in string_tfs){
  row = which(rownames(delta_matrix_24v0)==toupper(tf))
  rows = c(rows,row)
}
out_crop=delta_matrix_24v0[rows,]




out_order_24v0=order(out_crop[,3])

pdf(file='C:/Users/rhirsch/Documents/out_degree_log2_fc_brd4_24v0_error_bars_string_tfs.pdf',width=18)

barCenters <- barplot(out_crop[out_order_24v0,3],las=2,main='Rank order by log2 FC Brd4 out degree 24H vs 0H',ylim=c(-.2,.2))
segments(barCenters, out_crop[out_order_24v0,3] - out_crop[out_order_24v0,5], barCenters,out_crop[out_order_24v0,3] + out_crop[out_order_24v0,5], lwd = 1.5)
arrows(barCenters, out_crop[out_order_24v0,3] - out_crop[out_order_24v0,5], barCenters,out_crop[out_order_24v0,3] + out_crop[out_order_24v0,5], lwd = 1.5, angle = 90,code = 3, length = 0.05)

dev.off()

rows=c()
for(tf in string_tfs){
  row = which(rownames(delta_matrix_2v0)==toupper(tf))
  rows = c(rows,row)
}
out_crop_2v0=delta_matrix_2v0[rows,]

out_order_2v0=order(out_crop_2v0[,3])

pdf(file='C:/Users/rhirsch/Documents/out_degree_log2_fc_brd4_2v0_error_bars_string_tfs.pdf',width=18)

barCenters <- barplot(out_crop_2v0[out_order_2v0,3],las=2,main='Rank order by log2 FC Brd4 out degree 2H vs 0H',ylim=c(-.1,1.25))
segments(barCenters, out_crop_2v0[out_order_2v0,3] - out_crop_2v0[out_order_2v0,5], barCenters,out_crop_2v0[out_order_2v0,3] + out_crop_2v0[out_order_2v0,5], lwd = 1.5)
arrows(barCenters, out_crop_2v0[out_order_2v0,3] - out_crop_2v0[out_order_2v0,5], barCenters,out_crop_2v0[out_order_2v0,3] + out_crop_2v0[out_order_2v0,5], lwd = 1.5, angle = 90,code = 3, length = 0.05)

dev.off()


NHR_list = c('Thra','Thrb', 'Nr1d1', 'Rara', 'Rxra', 'Rxrb', 'Rarg', 'Esrra', 'Nr2f6','Nr2f2')
JF_list = c('Jun', 'Junb', 'Jund','Fosl1', 'Fosl2')
td_list= c('Tead1','Tead3')

pdf(file='C:/Users/rhirsch/Documents/stringdb_2H_24H_out_degree_scatter_jf_black_24_order_ef.pdf')
plot(as.numeric(out_crop_2v0[out_order_24v0,3]),as.numeric(out_crop[out_order_24v0,3]), col=ifelse(rownames(out_crop_2v0[out_order_24v0,]) %in% toupper(JF_list), "black", "grey"),
     pch=19,xlab = 'Log2 FC 2H vs 0H brd4 out degree', ylim=c(-.25,.25),xlim=c(.6,1.2),
     ylab = 'Log2 FC 24H vs 0H brd4 out degree',main='Comparison of Log2 24H v 0H Brd4 out degree and Log2 2H v 0H Brd4 out degree for string db TFs' )
dev.off()

pdf(file='C:/Users/rhirsch/Documents/stringdb_2H_24H_out_degree_scatter_nhr_red_24_order_ef.pdf')
plot(as.numeric(out_crop_2v0[out_order_24v0,3]),as.numeric(out_crop[out_order_24v0,3]), col=ifelse(rownames(out_crop_2v0[out_order_24v0,]) %in% toupper(NHR_list), "red", "grey"),
     ylim=c(-.25,.25),xlim=c(.6,1.2),pch=19,xlab = 'Log2 FC 2H vs 0H brd4 out degree',
     ylab = 'Log2 FC 24H vs 0H brd4 out degree',main='Comparison of Log2 24H v 0H Brd4 out degree and Log2 2H v 0H Brd4 out degree for string db TFs' )
dev.off()

pdf(file='C:/Users/rhirsch/Documents/stringdb_2H_24H_out_degree_scatter_td_blue_24_order_ef.pdf')
plot(as.numeric(out_crop_2v0[out_order_24v0,3]),as.numeric(out_crop[out_order_24v0,3]), col=ifelse(rownames(out_crop_2v0[out_order_24v0,]) %in% toupper(td_list), "blue", "grey"),
     ylim=c(-.25,.25),xlim=c(.6,1.2),pch=19,xlab = 'Log2 FC 2H vs 0H brd4 out degree',
     ylab = 'Log2 FC 24H vs 0H brd4 out degree',main='Comparison of Log2 24H v 0H Brd4 out degree and Log2 2H v 0H Brd4 out degree for string db TFs' )
dev.off()


rows=c()
for(tf in string_tfs){
  row = which(rownames(delta_matrix_24v2)==toupper(tf))
  rows = c(rows,row)
}
out_crop_24v2=delta_matrix_24v2[rows,]

out_order_24v2=order(out_crop_24v2[,3])

pdf(file='C:/Users/rhirsch/Documents/out_degree_log2_fc_brd4_24v2_error_bars_string_tfs.pdf',width=18)

barCenters <- barplot(out_crop_24v2[out_order_24v2,3],las=2,main='Rank order by log2 FC Brd4 out degree 24H vs 2H',ylim=c(-1.5,.1),names=rownames(out_crop_24v2[out_order_24v2,]))
segments(barCenters, out_crop_24v2[out_order_24v2,3] - out_crop_24v2[out_order_24v2,5], barCenters,out_crop_24v2[out_order_24v2,3] + out_crop_24v2[out_order_24v2,5], lwd = 1.5)
arrows(barCenters, out_crop_24v2[out_order_24v2,3] - out_crop_24v2[out_order_24v2,5], barCenters,out_crop_24v2[out_order_24v2,3] + out_crop_24v2[out_order_24v2,5], 
                lwd = 1.5, angle = 90,code = 3, length = 0.05)

dev.off()
