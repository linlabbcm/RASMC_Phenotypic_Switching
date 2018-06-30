#setwd('~/Dropbox/rasmc/')
install.packages('plotrix')

library(plotrix)
#==================================================================
#========================HELPER FUNCTIONS==========================
#==================================================================

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
        }

#==================================================================
#========================DATA INPUT================================
#==================================================================

top_table_24 = read.delim('/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_24H_top_table_1000.txt')
bottom_table_24 = read.delim('/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_24H_bottom_table_1000.txt')

top_table_2 = read.delim('/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_2H_top_table_1000.txt')
bottom_table_2 = read.delim('/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_2H_bottom_table_1000.txt')



#==================================================================
#=========COMPARING TF ENRICHMENT AT TOP VS BOTTOM 1000 GENES======
#==================================================================



plot_enrichment <- function(top_table,bottom_table,nSample=1000){
	tf_list = as.character(colnames(top_table)[7:ncol(top_table)])
	
	enrichment_vector = c()
	ue_vector = c()
	le_vector = c()

	top_length = sum(top_table$SUBPEAK_LENGTH)
	bottom_length = sum(bottom_table$SUBPEAK_LENGTH)
	
	for(tf in tf_list){
		print(tf)
		tf_col = which(colnames(top_table) == tf)
		tf_enrich_vector = c()
		for(n in 1:nSample){
			sample_rows_top = sample(1:nrow(top_table),nrow(top_table),TRUE)
			top_density = sum(top_table[sample_rows_top,tf_col])/sum(top_table[sample_rows_top,5])
			#top_density = sum(top_table[,tf_col])/sum(top_table[,5])
		
			sample_rows_bottom = sample(1:nrow(bottom_table),nrow(bottom_table),TRUE)
			bottom_density = sum(bottom_table[sample_rows_bottom,tf_col])/sum(bottom_table[sample_rows_bottom,5])
			tf_enrich_vector = c(tf_enrich_vector,-1 * log2(top_density/bottom_density))
			
		}
	
		enrichment_vector = c(enrichment_vector,median(tf_enrich_vector))
		ue_vector = c(ue_vector,quantile(tf_enrich_vector,0.975))
		le_vector = c(le_vector,quantile(tf_enrich_vector,0.025))
		
	}
	
	enrich_order = order(enrichment_vector)
	enrich_bar = barplot(enrichment_vector[enrich_order],names = tf_list[enrich_order],las=2,ylim=c(-0.7,0.7))
	error.bar(enrich_bar,enrichment_vector[enrich_order],enrichment_vector[enrich_order]-le_vector[enrich_order],ue_vector[enrich_order]-enrichment_vector[enrich_order])
}


pdf(file='/storage/cylin/grail/projects/rasmc_all/TF_Network/brd4_log2_24_v_0_ordered_top_1000_bar.pdf',width = 15,height = 5)
plot_enrichment(top_table_24,bottom_table_24,10000)
dev.off()


pdf(file='/storage/cylin/grail/projects/rasmc_all/TF_Network/brd4_log2_2_v_0_ordered_top_1000_bar.pdf',width = 15,height = 5)
plot_enrichment(top_table_2,bottom_table_2,10000)
dev.off()
