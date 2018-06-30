library(stats)

jf_table=read.delim('C:/Users/rhirsch/Documents/srf_analysis/subpeak_to_region_motif_density_table_jf.txt')
carg_table = read.delim('C:/Users/rhirsch/Documents/srf_analysis/subpeak_to_region_motif_density_carg_box_table.txt')

jf_cut_rows = which(jf_table$MOTIFS.KB >= .3)

carg_rows= which(carg_table$MOTIF_COUNT>0)
jf_rows = which(jf_table$MOTIF_COUNT>0)

carg_pos_table=carg_table[carg_rows,]
jf_pos_table=jf_table[jf_rows,]

#brd4_cut_rows = which(jf_table$LOG2FC_2v0 >= .585)

#cut_rows = intersect(jf_cut_rows,brd4_cut_rows)

crop_table = jf_table[jf_cut_rows,]


rows=c()
for(i in 1:length(crop_table[,1])){
  row=which(carg_table[,1]==crop_table[i,1])
  rows=c(rows,row)
  if(i%%1000 == 0){
  print(i)
  }
}

carg_crop_table=carg_table[rows,]


zero_srf_rows=which(carg_crop_table[,7]==0)
srf_rows=which(carg_crop_table[,7]!=0)



high_srf=carg_crop_table[srf_rows,]
low_srf=carg_crop_table[zero_srf_rows,]

pdf(file='C:/Users/rhirsch/Documents/srf_analysis/srf_analyis_brd4_boxplots.pdf',width=18)
boxplot(jf_table$BRD4_UNSTIM,jf_table$BRD4_2H,jf_table$BRD4_24H,ylim=c(0,1),
        names=c('0H','2H','24'),outline=FALSE,main='Brd4 signal across timepoints at all subpeaks',ylab='Brd4 signal')
boxplot(carg_pos_table$BRD4_UNSTIM,jf_pos_table$BRD4_UNSTIM,carg_pos_table$BRD4_2H,jf_pos_table$BRD4_2H,carg_pos_table$BRD4_24H,jf_pos_table$BRD4_24H,ylim=c(0,1),
        outline=FALSE, main='Brd4 signal across high srf motif subpeaks and jf motif subpeaks',names=c('0H srf','0H jf','2H srf', '2H jf', '24H srf','24H jf'))
boxplot(high_srf$BRD4_UNSTIM,low_srf$BRD4_UNSTIM,high_srf$BRD4_2H,low_srf$BRD4_2H,high_srf$BRD4_24H,low_srf$BRD4_24H,ylim=c(0,1),
        names=c('0H high srf','0H low srf','2H high srf','2H low srf','24H high srf','24H low srf'),outline=FALSE,
        main='Brd4 signal across timepoints at JF motif enriched subpeaks with high and low Srf motif enrichment',ylab='Brd4 signal')
boxplot(high_srf$LOG2FC_2v0,low_srf$LOG2FC_2v0,high_srf$LOG2FC_24v0,low_srf$LOG2FC_24v0,ylim=c(-2,2),
        names=c('Log2v0 high srf','Log2v0 low srf','Log24v0 high srf','Log24v0 low srf'),outline=FALSE,
        main='Log2 FC Brd4 signal across timepoints at JF motif enriched subpeaks with high and low Srf motif enrichment',ylab='Log2 FC Brd4 signal')
dev.off()

t.test(high_srf$BRD4_UNSTIM,low_srf$BRD4_UNSTIM)
t.test(high_srf$BRD4_2H,low_srf$BRD4_2H)
t.test(high_srf$BRD4_24H,low_srf$BRD4_24H)

t.test(high_srf$LOG2FC_2v0,low_srf$LOG2FC_2v0)
t.test(high_srf$LOG2FC_24v0,low_srf$LOG2FC_24v0)












# 
# nBins=20
# nIter = 1000
# v1=carg_crop_table$LOG2FC_2v0
# v1Order = order(v1)
# v2=carg_crop_table$MOTIFS.KB
# 
# binSize = length(v1)/nBins
# binMatrix = matrix(ncol = nBins,nrow=binSize)
# i = 1
# for(i in 1:nBins){
#   start = 1 + (i-1)*binSize
#   stop = i*binSize
#   
#   binMatrix[,i] = as.numeric(v2[v1Order[start:stop]])
#   
#   
# }
# 
# crop_bin_mean_vector=apply(binMatrix,2,mean)
# 
# meanMatrix = matrix(ncol = nBins,nrow=nIter)
# 
# for(i in 1:nIter){
#   for(j in 1:nBins){
#     meanMatrix[i,j] = mean(sample(binMatrix[,j],binSize,replace=TRUE))
#   }
# }
# 
# 
# meanVector = apply(meanMatrix,2,mean)
# upperError = apply(meanMatrix,2,quantile,probs=0.975) - apply(meanMatrix,2,mean)
# lowerError = apply(meanMatrix,2,mean) - apply(meanMatrix,2,quantile,probs=0.025)
# 
# 
# 
# crop_mean_mat_vector = apply(meanMatrix,2,mean)
# 
# 
# 
# 
# #######################################################
# ###########################ALLLLL######################
# #######################################################
# 
# nBins=680
# 
# v1=carg_table$LOG2FC_2v0
# v1Order = order(v1)
# v2=carg_table$MOTIFS.KB
# 
# binSize = length(v1)/nBins
# binMatrix = matrix(ncol = nBins,nrow=binSize)
# i = 1
# for(i in 1:nBins){
#   start = 1 + (i-1)*binSize
#   stop = i*binSize
#   
#   binMatrix[,i] = as.numeric(v2[v1Order[start:stop]])
#   
#   
# }
# 
# carg_bin_mean_vector=apply(binMatrix,2,mean)
#   
# meanMatrix = matrix(ncol = nBins,nrow=nIter)
# 
# for(i in 1:nIter){
#   for(j in 1:nBins){
#     meanMatrix[i,j] = mean(sample(binMatrix[,j],binSize,replace=TRUE))
#   }
# }
# 
# 
# meanVector = apply(meanMatrix,2,mean)
# upperError = apply(meanMatrix,2,quantile,probs=0.975) - apply(meanMatrix,2,mean)
# lowerError = apply(meanMatrix,2,mean) - apply(meanMatrix,2,quantile,probs=0.025)
# 
# 
# carg_mean_mat_vector = apply(meanMatrix,2,mean)
# 
# 
# 
# pdf(file='C:/Users/rhirsch/Documents/srf_analysis/carg_box_boxplots_under_JF_and_ALL_subpeaks.pdf')
# boxplot(crop_mean_mat_vector,carg_mean_mat_vector,names=c('Jun/Fos subpeaks','All subpeaks'),ylab='Mean CArG box motifs/kbp',main='Mean CArG box motifs/kbp across subpeaks')
# dev.off()
# 
# 
