options(scipen=999)

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
  colorSpectrum <- colorRampPalette(c("blue","grey","grey","red"))(100)
  
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
#pdf(file = "/home/rahirsch/Desktop/test_pdf.pdf")
#v1 = jitter(1:1000,amount = 100) #makes a jittery vector that goes from 1 to 1000
#v2 = jitter((1:1000)^2,amount = 100) #makes a jittery exponential vector

#compareVectors(v1,v2,'v1','v2','example')
#dev.off()

region_table = read.delim('/storage/cylin/grail/projects/rasmc_all/motif_density/subpeak_coords_table.txt',header=FALSE)
dim(region_table)
signal_table = read.delim('/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_subpeak_Brd4_signal_table.txt',header=TRUE)
dim(signal_table)
density_table = read.delim('/storage/cylin/grail/projects/rasmc_all/motif_density/CArG_box_seq_density_from_fasta_full_length.txt',header=TRUE)
dim(density_table)
density_table[1:5,]
region_map = read.delim('/storage/cylin/grail/projects/rasmc_all/TF_Network/regions_subpeak_order.txt',header=FALSE)
dim(region_map)


#Load in tables
mean_unstim=c()
for(i in 1:length(signal_table[,1])){
  mean = sum(signal_table[i,9],signal_table[i,10])/2
  mean_unstim=c(mean_unstim,mean)
}

mean_2h=c()
for(i in 1:length(signal_table[,1])){
  mean = sum(signal_table[i,11],signal_table[i,12])/2
  mean_2h=c(mean_2h,mean)
}

mean_24h=c()
for(i in 1:length(signal_table[,1])){
  mean = sum(signal_table[i,13],signal_table[i,14])/2
  mean_24h=c(mean_24h,mean)
}



compare_table = cbind(region_table[,1:4],region_map[,1],density_table[,1]*1000,round(density_table[,1]*density_table[3]),density_table[,2:3],mean_unstim,mean_2h,mean_24h,log2(mean_2h/mean_unstim),log2(mean_24h/mean_unstim))

colnames(compare_table)=c('SUBPEAK_ID','CHROM','START','STOP','REGION_ID','MOTIFS/KB','MOTIF_COUNT','POSITIONS','SUBPEAK_LENGTH','BRD4_UNSTIM','BRD4_2H','BRD4_24H','LOG2FC_2v0','LOG2FC_24v0')

compare_table[1:5,]

write.table(compare_table,'/storage/cylin/grail/projects/rasmc_all/motif_density/subpeak_to_region_motif_density_carg_box_table.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

#example
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/compare_vectors_log2brd4_24v0_carg_box_motif_density.pdf")
v2 = compare_table$`MOTIFS/KB`
v1 = compare_table$LOG2FC_24v0

compareVectors(v1,v2,'brd4 24_v_0','CArG Box motif density','Log2 Fold Change in 24H Brd4 versus motif density')
dev.off()



#example
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/compare_vectors_log2brd4_2v0_carg_box_motif_density.pdf")
v2 = compare_table$`MOTIFS/KB`
v1 = compare_table$LOG2FC_2v0

compareVectors(v1,v2,'brd4 2_v_0','CArG Box motif density','Log2 Fold Change in 2H Brd4 versus motif density')
dev.off()


#example
pdf(file = "/storage/cylin/grail/projects/rasmc_all/figures/compare_vectors_log2brd4_2v24_carg_box_motif_density.pdf")
v2 = compare_table$`MOTIFS/KB`
v1 = log2(mean_24h/mean_2h)

compareVectors(v1,v2,'brd4 2_v_24','CArG Box motif density','Log2 Fold Change in 2H v 24H Brd4 versus motif density')
dev.off()


