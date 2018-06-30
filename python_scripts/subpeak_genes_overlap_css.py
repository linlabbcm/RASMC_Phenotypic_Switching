import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils

subpeak_genes_path = '/storage/cylin/grail/projects/rasmc_all/crc/rasmc_h3k27ac_0_tss/subpeak_genes_regions.txt'

subpeak_signal_table = '/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_subpeak_Brd4_signal_table.txt'

subpeak_genes = utils.parseTable(subpeak_genes_path, '\t')
subpeak_signal = utils.parseTable(subpeak_signal_table,'\t')

overlap_table = []

for line in subpeak_signal[1:]:
	gene_overlaps=[]
	chrom = line[1].split('(')[0]
	tss_1 = line[1].split(':')[1]
	tss_sig = int(tss_1.split('-')[0])
	stop_sig = int(tss_1.split('-')[1])
	for gene in subpeak_genes:
		chr = gene[1]
   		start = int(gene[2])
		stop = int(gene[3])
		gName = gene[0]
		if chrom == chr:
			if start>=tss_sig and tss_sig<=stop:
				gene_overlaps.append(gName)
			if stop<=stop_sig and stop>=stop_sig:
				gene_overlaps.append(gName)
	gene_string =  ",".join(gene_overlaps)
	new_line = [gene_string,chrom,tss_sig,stop_sig,line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12]]
	overlap_table.append(new_line)

utils.unParseTable(overlap_table,'/storage/cylin/grail/projects/rasmc_all/TF_Network/subpeak_genes_signal_table_css.txt','\t')
