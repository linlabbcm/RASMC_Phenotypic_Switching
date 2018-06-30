import sys
sys.path.append('/storage/cylin/bin/pipeline/')

import pipeline_dfci
import utils
import string

from collections import defaultdict

genes_list_path = '/storage/cylin/grail/projects/rasmc_all/TF_Network/crc_genes_list.txt'
regions_24H_path = '/storage/cylin/grail/projects/rasmc_all/TF_Network/170206_waterfall_genes_log2_24_v_0_ordered_regions.txt'

regions_2H_path = '/storage/cylin/grail/projects/rasmc_all/TF_Network/170206_waterfall_genes_log2_2_v_0_ordered_regions.txt'

unique_regions_path = '/storage/cylin/grail/projects/rasmc_all/TF_Network/unique_regions_subpeaks_lengths.txt'

tf_targets_path = '/storage/cylin/grail/projects/rasmc_all/TF_Network/tf_targets_list_table.txt'

genes_list = utils.parseTable(genes_list_path,'\t')
regions_24H = utils.parseTable(regions_24H_path,'\t')
regions_2H = utils.parseTable(regions_2H_path,'\t')
unique_regions = utils.parseTable(unique_regions_path,'\t') 
tf_targets = utils.parseTable(tf_targets_path,'\t')

genesDict = defaultdict(list)

print(regions_24H[1])

for line in regions_24H[1:]:
	genesDict[line[1].upper()]=[line[4]]

targetDict = defaultdict(list)

for target in tf_targets[1:]:
	targetDict[target[0]]=[target[1]]



print(genesDict['TEAD3'])
print(targetDict['TEAD3'])
overlap_matrix = []

header = ['PEAK_ID','CHROM','START','STOP','SUBPEAK_LENGTH','GENE_LIST']

for gene in genes_list:
	header.append(gene[0])

overlap_matrix.append(header)

target_matrix = []
target_header=['target_gene','regions']
for gene in genes_list:
	target_header.append(gene[0])
target_matrix.append(target_header)

for line in regions_24H[1:]:
	target_gene = line[1].upper()
	overlap_list=[]
	for gene in genes_list:
		if any(target_gene in target for target in targetDict[gene[0]]):
			overlap_list.append(1)
		else:
			overlap_list.append(0)
#        if any('1' in str(overlap) for overlap in overlap_list):
        new_line=[target_gene,line[4]]
        new_line.append(overlap_list)
        target_matrix.append(new_line)

print(target_matrix[0:10])

derp = 1

 for line in unique_regions:
 	print(derp)
 	region = line[0]
 	sp_len = line[1]
 	chrom = ''
 	start = ''
 	stop = ''
 	gene_list=[]
 	for foo in target_matrix:
 		if any(region in r for r in foo[1].split(',')):
 			gene_list.append(foo[0])	
 			new_line = [region,chrom,start,stop,sp_len,','.join(gene_list)]
 			new_line.extend(foo[2:110])
 			overlap_matrix.append(new_line)
 	derp=derp+1
 print(overlap_matrix[0:10])
 print(len(overlap_matrix))

utils.unParseTable(overlap_matrix,'/storage/cylin/grail/projects/rasmc_all/TF_Network/regions_genes_overlap_matrix.txt','\t')

def cut_1000(order_table,regions_table,geneDict, outpath_top='',outpath_bottom=''):

	top_table = []
	top_cut_regions = []
	top_table.append(regions_table[0])
	for line in order_table[1:1001]:
		top_cut_regions.append(line[4].split(','))
	print(top_cut_regions[1:10])
	for line in regions_table:
		region = line[0]
				#print(region)
		if any(region in s for s in top_cut_regions):
			top_table.append(line)


	utils.unParseTable(top_table,outpath_top,'\t')

	bottom_table=[]
	bottom_cut_regions = []
	bottom_table.append(regions_table[0])
	for line in order_table[len(order_table)-1000:len(order_table)]:
		bottom_cut_regions.append(line[4].split(','))

	print(bottom_cut_regions[1:10])
	for line in regions_table:
		region = line[0]
		if any(region in s for s in bottom_cut_regions):
			bottom_table.append(line)

	utils.unParseTable(bottom_table,outpath_bottom,'\t')

overlap_matrix = utils.parseTable('/storage/cylin/grail/projects/rasmc_all/TF_Network/regions_genes_overlap_matrix.txt','\t')

print(overlap_matrix[0:5])

olm_new = []
olm_new.append(overlap_matrix[0])
for line in overlap_matrix[1:]:
	r = line[0]
	c = line[1]
	start = line[2]
	stop = line[3]
	sp_len =line[4]
	genes = line[5]
	counts = ((line[6].split('[')[1]).split(']')[0]).split(', ')
	new_line = [r,c,start,stop,sp_len,genes]
	for count in counts:
		new_line.append(int(count))
	olm_new.append(new_line)

overlap_matrix = olm_new

print(overlap_matrix[0:5])
cut_1000(regions_24H,overlap_matrix,genesDict,'/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_24H_top_table_1000.txt','/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_24H_bottom_table_1000.txt')


cut_1000(regions_2H,overlap_matrix,genesDict,'/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_2H_top_table_1000.txt','/storage/cylin/grail/projects/rasmc_all/TF_Network/enrichment_2H_bottom_table_1000.txt')
