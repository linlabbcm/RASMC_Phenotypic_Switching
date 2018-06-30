import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils
from collections import defaultdict

expr_path = '/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_all_fpkm_means.txt'

active_genes_path='/storage/cylin/grail/projects/rasmc_all/tables/activeListTableGENES.txt'
rat_path='/storage/cylin/bin/pipeline/crc/annotation/TFlist_NMid_rn6.txt'

rat_tfs = utils.parseTable(rat_path,'\t')
exp_table= utils.parseTable(expr_path, '\t')
active_table=utils.parseTable(active_genes_path,'\t')

exp_dict=defaultdict(list)

for i in exp_table[1:]:
	print(i)
	gene = i[0]
#	print(gene)
	exp0 = float(i[1])
#	print(i[1])
	exp2 = float(i[2])
#	print(i[2])
	exp24=float(i[4])
#	print(i[4])
	exp_dict[gene]=[exp0,exp2,exp24]


ag_cut_list=[]

for gene in active_table:
	print(gene)
	foo=exp_dict[gene[0]]
	if len(foo) == 3:
		if foo[0] >= 10 or foo[1] >=10 or foo[2] >= 10:
			ag_cut_list.append(gene)


print(len(ag_cut_list))

tfs=[]

for gene in ag_cut_list:
	count=0
	name=gene[0].upper()
	for tf in rat_tfs:
		if name == tf[1]:
			count=count+1
	if count > 0:
		tfs.append(gene[0])

print(len(tfs))