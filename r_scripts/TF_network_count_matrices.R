#load in tables
edgeTable = read.delim('/storage/cylin/grail/projects/rasmc_all/tables/unique_filtered_edge_exp_closest_genes_filter.txt',header=FALSE)

TF_list = edgeTable[,1]
target_list = edgeTable[,2]

TF_list = unique(TF_list)
target_list = unique(target_list)

o=order(target_list)

target_list=target_list[o]

countMatrix = matrix(nrow=length(TF_list),ncol=length(target_list))
rownames(countMatrix)=TF_list
colnames(countMatrix)=target_list

t_rows = c()
count_rows=c()

for(i in 1:length(TF_list)){
  print(i)
  t_rows = which(edgeTable[,1]==as.character(TF_list[i]))
  for(j in 1:length(target_list)){

    print(as.character(target_list[j]))
    count_rows = which(edgeTable[t_rows,2] == as.character(target_list[j]))
    count = length(count_rows)
    
    
    countMatrix[i,j] = count  
    
  }
}

write.table(countMatrix,file='/storage/cylin/grail/projects/rasmc_all/TF_Network/tf_targets_count_matrix.txt',quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)


binaryCountMatrix = matrix(nrow=length(TF_list),ncol=length(target_list))
rownames(binaryCountMatrix)=TF_list
colnames(binaryCountMatrix)=target_list

for(i in 1:nrow(countMatrix)){
  for(j in 1:ncol(countMatrix)){
    if(countMatrix[i,j] > 0){
      binaryCountMatrix[i,j] = 1
    }
    else{
      binaryCountMatrix[i,j] = 0
    }
  }
}

write.table(binaryCountMatrix,file='/storage/cylin/grail/projects/rasmc_all/TF_Network/tf_targets_binary_count_matrix.txt',quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)


