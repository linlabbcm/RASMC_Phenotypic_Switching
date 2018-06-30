table = read.delim('C:/Users/rhirsch/Desktop/rasmc_enhancer_promoter/brd4_tables/subpeak_to_region_motif_density_table.txt', header=TRUE)

scatter_table = cbind(table[,1:5],table[,10:14],log2(table[,12]/table[,11]))
colnames(scatter_table)=c("SUBPEAK_ID","CHROM","START","STOP","REGION_ID","BRD4_UNSTIM","BRD4_2H","BRD4_24H","LOG2FC_2v0","LOG2FC_24v0",'LOG2FC_24v2')


plot(scatter_table[,9],scatter_table[,11], xlim=c(-3,3),ylim=c(-3,3),xlab="Log2FC 2H versus Unstim", ylab="Log2FC 24H versus 2H", pch=1, cex=1)
abline(h=0, v=0, lwd=2)

#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================

library(ggplot2)


vector1 = scatter_table[,9]
vector2 = scatter_table[,11]

#==================================================================
#=========================DENSITY PLOTS============================
#==================================================================

#for some reason ggplot doesn't like looping

print('Plotting contour scatter plot for RASMC subpeak regions')
pdf(file='C:/Users/rhirsch/Desktop/rasmc_enhancer_promoter/brd4_tables/subpeak_to_region_scatter_plot_2v0_vs24v2.pdf',width =6.5,height =5)
dataset = structure(list(x= vector1,y=vector2),class = "data.frame")


ggplot(dataset, aes(x, y)) + 
  geom_point(data = dataset,cex=0.5,pch=16,col=rgb(0.5,0.5,0.5,0.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=50, geom="polygon") +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
  geom_density2d(colour=rgb(0.5,0.5,0.5,.8)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))





dev.off()



vector1 = scatter_table[,9]
vector2 = scatter_table[,10]


print('Plotting contour scatter plot for RASMC subpeak regions')
pdf(file='C:/Users/rhirsch/Desktop/rasmc_enhancer_promoter/brd4_tables/subpeak_to_region_scatter_plot_2v0_vs24v0.pdf',width =6.5,height =5)
dataset = structure(list(x= vector1,y=vector2),class = "data.frame")


ggplot(dataset, aes(x, y)) + 
  geom_point(data = dataset,cex=0.5,pch=16,col=rgb(0.5,0.5,0.5,0.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=2, bins=50, geom="polygon") +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0.00, 0.5), guide = FALSE) +
  geom_density2d(colour=rgb(0.5,0.5,0.5,.8)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))





dev.off()    

