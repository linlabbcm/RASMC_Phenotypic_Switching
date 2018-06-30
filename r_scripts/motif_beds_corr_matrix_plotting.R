options(scipen=999)
library(GenomicFeatures)



args[1] = wkdir
args[2] = tf_path
args[3] = analysis_name
args[4] = window

tf_list_table=read.delim(tf_path,header=F,sep='\t')
tf_list=c(as.character(window),tf_list_table[,1])



all.peaks.files <- list.files(path = paste(wkdir,'tmp/',sep='') , pattern = '*.bed$',  full.names = T)


tf_files=c(1)
for(j in 1:length(tf_list)){
  for(i in 1:length(all.peaks.files)){
    fn = unlist(strsplit(toString(all.peaks.files[i]),'/tmp/'))
    foo = fn[2]
    name = unlist(strsplit(toString(foo),'_'))
    if(tf_list[j]==name[1]){
      tf_files=c(tf_files,i)
    }
  }
}  

tf.peaks.files <- all.peaks.files[tf_files]

ls.gr.enhancers <- lapply(tf.peaks.files, function(file.name) {
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
  print(i)
}

names=c()
for(i in 1:length(all.peaks.files)){
  file_name = unlist(strsplit(toString(all.peaks.files[i]),'/tmp/'))
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

corrMatrix = as.matrix(overlapMatrix[,7:length(overlapMatrix)])
rownames(corrMatrix) = row_names

write.table(corrMatrix, file = paste(wkdir,'tables/','corr_matrix_',window,'.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

keep_rows = c()
for(i in 1:length(corrMatrix[,1])){
  count = 0
  for(j in 1:length(corrMatrix[1,])){
    if(corrMatrix[i,j] == 1){
      count = count + 1
    }
  }
  if(count > 1){
    keep_rows = c(keep_rows,i)
  }
}

length(keep_rows)

fcorrMatrix = corrMatrix[keep_rows,]


#for 50bp

tf_matrix_path = paste(wkdir,'tables/','corr_matrix_',window,'.txt',sep='')

tf_matrix = read.delim(tf_matrix_path)


#for euclidean distance
dist_matrix = dist(t(tf_matrix))
tf_clust = hclust(dist_matrix)


#made a matrix where rows and columns are TF names
tf_cor_matrix = matrix(nrow=ncol(tf_matrix),ncol=ncol(tf_matrix))
rownames(tf_cor_matrix) = colnames(tf_matrix)
colnames(tf_cor_matrix) = colnames(tf_matrix)


for(i in 1:ncol(tf_matrix)){
  
  for(j in 1:ncol(tf_matrix)){
    
    #set distance to 1- pearson correlation
    tf_cor_matrix[i,j] = as.numeric(1-cor(tf_matrix[,i],tf_matrix[,j]))
    
  }
  
  
}

#made a matrix where rows and columns are TF names
tf_c_matrix = matrix(nrow=ncol(tf_matrix),ncol=ncol(tf_matrix))
rownames(tf_c_matrix) = colnames(tf_matrix)
colnames(tf_c_matrix) = colnames(tf_matrix)

for(i in 1:ncol(tf_matrix)){
  
  for(j in 1:ncol(tf_matrix)){
    
    #pearson correlation
    tf_c_matrix[i,j] = as.numeric(cor(tf_matrix[,i],tf_matrix[,j]))
    
  }
  
  
}


#tf_clust = hclust(as.dist(tf_cor_matrix))
#plot(tf_clust)


#==============================================================================
#==============================================================================
#==============================================================================


#===================================================================
#====================MAKING ENHANCER MATRIX=========================
#===================================================================

enhancerTable = tf_c_matrix




enhancerMatrix = as.matrix(enhancerTable)
#rownames(enhancerMatrix) = enhancerTable[,1]
#===================================================================
#==================SAMPLE DISTANCE MATRICIES========================
#===================================================================

#distance by samples
sampleDist = as.dist(1-cor(enhancerMatrix))
sampleHC = hclust(sampleDist)


#===================================================================
#=====================PLOTTING SAMPLE TREE==========================
#===================================================================


#PLOTTING THE SAMPLE TREE

treePDFFile = paste(wkdir,'figures/',analysisName,"_treePlot.pdf",sep='')
pdf(treePDFFile,width=8,height=8)


plot(sampleHC,xlab='',main = paste(genome,'_',analysisName,sep=''))
dev.off()



#===================================================================
#=======================SAMPLE CLUSTER ORDERING=====================
#===================================================================


sampleOrder = sampleHC$order



#===================================================================
#=================MAKING SAMPLE PAIRWISE HEATMAP====================
#===================================================================



#distance by samples
sampleSimMatrix = cor(enhancerMatrix)


#Set the color spectrum
colorSpectrum <- colorRampPalette(c("blue","white","red"))(100)

#setting a color data range
minValue <- -1
maxValue <- 1
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(0, color_cuts,1)

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)


clusterSampleFile = paste(wkdir,'figures/',analysisName,"_clusterSamples_range_dist.pdf",sep='')

pdf(file = clusterSampleFile,width = 10,height =10)

layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))
par(mar=c(8,9,5,9))
plot(as.dendrogram(sampleHC),ylab='Distance')
par(mar=c(4,2,2,2))

plot(as.dendrogram(sampleHC),horiz=TRUE,xlab='Distance',leaflab='none')
par(mar=c(6,4,4,2))

image(1:ncol(sampleSimMatrix),1:nrow(sampleSimMatrix),t(sampleSimMatrix[sampleOrder,sampleOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')

par(mar=c(6,5,4,2))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Similarity")
dev.off()


