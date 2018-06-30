#161027_rasmc_edge_weights.R

#setwd('~/Dropbox/rasmc/')

#add cutoffs for expression
#=========================================================
#========================HELPER FUNCTIONS=================
#=========================================================

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
        }


print('Loading edge table')
edgeTable = read.delim('/storage/cylin/grail/projects/rasmc_all/crc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_EDGE_TABLE.txt')

print('Loading Brd4 signal table')
brd4Table = read.delim('/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_EDGE_TABLE_Brd4_signal_table.txt')

print('Loading expression table')
expTable = read.delim('/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_exprs_fpkm_means.txt')


print('Making out matrix')
brd4_unstim = (brd4Table[,9] + brd4Table[,10])/2
brd4_2h = (brd4Table[,11] + brd4Table[,12])/2
brd4_24h = (brd4Table[,13] + brd4Table[,14])/2

out_matrix = as.matrix(cbind(brd4_unstim,brd4_2h,brd4_24h))
rownames(out_matrix) = as.character(edgeTable[,1])

plot(ecdf(out_matrix[,1]),log='x',xlim = c(0.05,5))
lines(ecdf(out_matrix[,2]),col='blue')
lines(ecdf(out_matrix[,3]),col='red')


#median normalize
#boxplot(out_matrix[,1],out_matrix[,2],out_matrix[,3])

median_out_matrix = as.matrix(out_matrix)

median_out_matrix[,1] =median_out_matrix[,1]/median(median_out_matrix[,1])
median_out_matrix[,2] =median_out_matrix[,2]/median(median_out_matrix[,2])
median_out_matrix[,3] =median_out_matrix[,3]/median(median_out_matrix[,3])

#doing the 0 to 24 h edge comparison
plot(ecdf(median_out_matrix[,1]),log='x',xlim = c(0.1,10),lwd=2)
lines(ecdf(median_out_matrix[,2]),col='blue',lwd=2)
lines(ecdf(median_out_matrix[,3]),col='red',lwd=2)


top_rows = which(apply(median_out_matrix,1,max) > 1)

filtered_out_matrix = median_out_matrix[top_rows,]

makeFoldMatrix <- function(m,expTable,nIter = 10){

	tf_names = unique(rownames(m))
	
	exp_tf_names = c()
	for(tf in tf_names){
		row = which(toupper(rownames(expTable))==tf)
		if(max(expTable[row,c(1,2,4)]) >= 10){
			exp_tf_names = c(exp_tf_names,tf)
			
		}
		
	}
	print(exp_tf_names)
	
	fold_matrix  = matrix(nrow = length(exp_tf_names),ncol = 6 )
	
	
	for(i in 1:length(exp_tf_names)){
		tf 	= exp_tf_names[i]
		tf_rows = which(rownames(m)==tf)
		print(tf)
		print(length(tf_rows))
		
		#for 2 hours
		mean_change = log2(m[tf_rows,2]/m[tf_rows,1])
		mean_change[which(mean_change <= -6)]<- -6 #constrain fold changes to reasonable numbers
		mean_change[which(mean_change >= 6)]<- 6	
		
		#now do resampling		
		mean_vector = c()
		
		for( n in 1:nIter){
			
			mean_vector = c(mean_vector,mean(sample(mean_change,length(mean_change),replace=TRUE)))
			
		}
		deets = quantile(mean_vector,probs = c(0.025,.5,0.975))
		fold_matrix[i,1:3] = deets
				
		#for the 24 hour time point
		mean_change = log2(m[tf_rows,3]/m[tf_rows,1])
		mean_change[which(mean_change <= -6)]<- -6 #constrain fold changes to reasonable numbers
		mean_change[which(mean_change >= 6)]<- 6	
		
		#now do resampling		
		mean_vector = c()
		
		for( n in 1:nIter){
			
			mean_vector = c(mean_vector,mean(sample(mean_change,length(mean_change),replace=TRUE)))
			
		}
		deets = quantile(mean_vector,probs = c(0.025,.5,0.975))
		fold_matrix[i,4:6] = deets		
	}
		

	rownames(fold_matrix) = exp_tf_names
	colnames(fold_matrix) = c('PDGF_2H_LOW','PDGF_2H_MEAN','PDGF_2H_HIGH','PDGF_24H_LOW','PDGF_24H_MEAN','PDGF_24H_HIGH')
	return(fold_matrix)
}




median_fold_matrix = makeFoldMatrix(median_out_matrix,expTable,1000)
filtered_fold_matrix = makeFoldMatrix(filtered_out_matrix,expTable,1000)

#writing the filtered_fold_matrix

order_24 = order(filtered_fold_matrix[,3])

write.table(filtered_fold_matrix[order_24,],file='/storage/cylin/grail/projects/rasmc_all/tables/RN6_RASMC_0_tss_TF_OUT_DEGREE_24_order.txt',quote=FALSE,sep='\t',row.names=TRUE,col.names=TRUE)


#============================

pdf(file='/storage/cylin/grail/projects/rasmc_all/figures/161130_rasmc_0_tss_tf_brd4_out_bar.pdf',width = 11,height=  8.5)
par(mfrow=c(2,1))

order_2 = order(filtered_fold_matrix[,2])

foo = barplot(filtered_fold_matrix[order_2,2],las=2,ylim = c(-0.2,0.4),main = 'PDGF 2 hour BRD4 out degree change',ylab='log2 change in BRD4 out degree')
error.bar(foo,filtered_fold_matrix[order_2,2],filtered_fold_matrix[order_2,3]-filtered_fold_matrix[order_2,2],filtered_fold_matrix[order_2,2]-filtered_fold_matrix[order_2,1])


order_24 = order(filtered_fold_matrix[,5])

foo = barplot(filtered_fold_matrix[order_24,5],las=2,ylim = c(-0.6,0.2),main = 'PDGF 24 hour BRD4 out degree change',ylab='log2 change in BRD4 out degree')
error.bar(foo,filtered_fold_matrix[order_24,5],filtered_fold_matrix[order_24,6]-filtered_fold_matrix[order_24,5],filtered_fold_matrix[order_24,5]-filtered_fold_matrix[order_24,4])

dev.off()





par(mfrow=c(2,1))

order_2 = order(filtered_fold_matrix[,1])
barplot(filtered_fold_matrix[order_2,1],las=2,main='2 hour filtered')
order_24 = order(filtered_fold_matrix[,2])
barplot(filtered_fold_matrix[order_24,2],las=2,main='24 hour filtered')



#=====================


par(mfrow=c(2,1))
order_2 = order(median_fold_matrix[,1])
barplot(median_fold_matrix[order_2,1],las=2,main='2 hour unfiltered')
order_24 = order(median_fold_matrix[,2])
barplot(median_fold_matrix[order_24,2],las=2,main ='24 hour unfiltered')



par(mfrow=c(2,1))
order_2 = order(median_fold_matrix[,1])
barplot(median_fold_matrix[order_2,1],las=2,main='2 hour unfiltered')
order_2 = order(filtered_fold_matrix[,1])
barplot(filtered_fold_matrix[order_2,1],las=2,main='2 hour filtered')

#fold_order_2 = order(out_fold_vector_2)


#fold_order_24 = order(out_fold_vector_24)
#plot(out_fold_vector[fold_order],type='h')

 