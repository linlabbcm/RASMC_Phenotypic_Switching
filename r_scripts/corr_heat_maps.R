library(GenomicFeatures)
library(gplots)


all.peaks.files <- list.files(path = '/storage/cylin/grail/projects/rasmc_all/crc/rasmc_h3k27ac_0_tss/motif_beds/sorted_beds' , pattern = '*.bed$',  full.names = T)
ls.gr.enhancers <- lapply(all.peaks.files, function(file.name) {
  print(file.name)
  df.ranges <- read.delim(file.name, sep = "\t", header = F)
  df.ranges <- df.ranges[ ,1:3]
  colnames(df.ranges) <- c("seqnames", "start", "end")
  gr.peaks <- as(df.ranges, "GRanges")
})
gr.enhancers <- unique(reduce(do.call(c, ls.gr.enhancers)))

overlapMatrix = as.data.frame(gr.enhancers)

for(i in 1:length(all.peaks.files)){
  qt.hits <- queryHits(findOverlaps(gr.enhancers, ls.gr.enhancers[[i]]))
  vector = rep(0,length(gr.enhancers))
  vector[qt.hits] = 1
  overlapMatrix = cbind(overlapMatrix,vector)
}

names=c()
for(i in 1:length(all.peaks.files)){
  file_name = unlist(strsplit(toString(all.peaks.files[i]),'/sorted_beds/'))
  foo = file_name[2]
  name = unlist(strsplit(toString(foo),'_'))
  names = c(names,name[1])
}


row_names = c()
for(i in 1:length(overlapMatrix[,1])){
  regID=paste(overlapMatrix[i,1],overlapMatrix[i,2],overlapMatrix[i,3], sep = " ",collapes=NULL)
  row_names = c(row_names,regID)
}

colnames(overlapMatrix) <- c(colnames(overlapMatrix[1:5]),names)
     
corrMatrix = as.matrix(overlapMatrix[,6:length(overlapMatrix)])
rownames(corrMatrix) = row_names



## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(corrMatrix), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(corrMatrix, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
## Plot heatmap 
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/test_corr_plot.pdf')
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 
dev.off()



