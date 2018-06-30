#!/usr/bin/python


import sys

sys.path.append('/storage/cylin/bin/pipeline/')
import utils
from collections import defaultdict
import numpy
#very quick script to add the BRD4 in degree to the out degree table

projectFolder = '/storage/cylin/grail/projects/rasmc_all/' #will have to set this appropriately for torta or your local


#this is the out degree table from the 161027_rasmc_edge_weights.R script
out_path = '%stables/RN6_RASMC_0_tss_TF_OUT_DEGREE.txt' % (projectFolder)

#this is the enhancer tf assignment from crc
enhancer_path = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_ENHANCER_TF_TABLE.txt' % (projectFolder)

#this is the signal table 
signal_path = '%ssignalTables/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE_signal_table.txt' % (projectFolder)

#this is the expression table
exp_path = '%scufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_exprs_fpkm_means.txt' % (projectFolder)

#this is the edge table path
edge_path = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_EDGE_TABLE.txt' % (projectFolder)


#load everybody in
out_table = utils.parseTable(out_path,'\t')

enhancer_table = utils.parseTable(enhancer_path,'\t')

signal_table = utils.parseTable(signal_path,'\t')

exp_table = utils.parseTable(exp_path,'\t')

edge_table = utils.parseTable(edge_path,'\t')

#first make an expression dictionary for each gene
#want each entry to be a list with the avg. exp at [0, 2+PDGF, 24+PDGF]

exp_dict = defaultdict(list)

for line in exp_table[1:]:
    exp_dict[line[0].upper()] = [float(line[1]),float(line[2]),float(line[4])]


#for each tf make a dict to hold AUC and then summarize in degree
#get this from the signal_table
brd4_auc = {'rasmc_unstim':defaultdict(float),
            'rasmc_2h':defaultdict(float),
            'rasmc_24h':defaultdict(float),
            }
for line in signal_table[1:]:
    
    gene_list = line[0].split(',')
    coords = [int(x) for x in line[1].split(':')[-1].split('-')]
    length = coords[1] - coords[0]
    for gene in gene_list:
        brd4_auc['rasmc_unstim'][gene]+=((float(line[8])+float(line[9]))/2)*length
        brd4_auc['rasmc_2h'][gene]+=((float(line[10])+float(line[11]))/2)*length
        brd4_auc['rasmc_24h'][gene]+=((float(line[12])+float(line[13]))/2)*length
        
        #sanity check for a gene w/ multiple SEs
        if gene == 'JUND':
            print('JUND')
            print(brd4_auc['rasmc_unstim']['JUND'])
            print(((float(line[8])+float(line[9]))/2)*length)

#sanity check for a pair of TFs that have the same SE
print('genes w/ same SE')
print(brd4_auc['rasmc_unstim']['NR1D1'])
print(brd4_auc['rasmc_unstim']['THRA'])


#now we have AUC for every single tf

#load in the out degree table

#set up the output table
pca_table_median = [['GENE','BRD4_DELTA_IN_2HR','BRD4_DELTA_IN_24HR','BRD4_DELTA_OUT_2HR','BRD4_DELTA_OUT_24HR','EXP_DELTA_IN_2HR','EXP_DELTA_IN_24HR','EXP_DELTA_OUT_2HR','EXP_DELTA_OUT_24HR']]

#get the median for each time point
#getting the median BRD4 signal across all TFs to median scale BRD4 IN Degree
brd4_unstim_median = numpy.median([brd4_auc['rasmc_unstim'][gene] for gene in brd4_auc['rasmc_unstim'].keys()])

brd4_2h_median = numpy.median([brd4_auc['rasmc_2h'][gene] for gene in brd4_auc['rasmc_2h'].keys()])

brd4_24h_median = numpy.median([brd4_auc['rasmc_24h'][gene] for gene in brd4_auc['rasmc_24h'].keys()])

print(brd4_unstim_median,brd4_2h_median,brd4_24h_median)


for line in out_table[1:]:
    gene = line[0]
    #just get the BRD4 out degree (nothing else)
    brd4_out_2 = float(line[2])
    brd4_out_24 = float(line[5])


    
    #ended up using just median scaled data
    brd4_in_2_median = numpy.log2((brd4_auc['rasmc_2h'][gene]/brd4_2h_median)/(brd4_auc['rasmc_unstim'][gene]/brd4_unstim_median))
    brd4_in_24_median = numpy.log2((brd4_auc['rasmc_24h'][gene]/brd4_24h_median)/(brd4_auc['rasmc_unstim'][gene]/brd4_unstim_median))

    #now getting the expression IN and OUT degree
    exp_in_2 = numpy.log2(exp_dict[gene][1]/exp_dict[gene][0])
    exp_in_24 = numpy.log2(exp_dict[gene][2]/exp_dict[gene][0])
    
    #for out degree we need to get all target genes
    #get all of the target edges from the edge table
    target_genes = utils.uniquify([line[1] for line in edge_table[1:] if line[0] == gene])

    #max exp values should be > 1
    # sad_genes = [target for target in target_genes if max(exp_dict[target]) < 1]
    # for sadness in sad_genes:
    #     print(exp_dict[sadness])
    # sys.exit()
    
    #now we can get the vector of log2 fold expression for each one
    #want the average exp change
    exp_out_2 = numpy.mean([numpy.log2(exp_dict[target][1]/exp_dict[target][0]) for target in target_genes])
    exp_out_24 = numpy.mean([numpy.log2(exp_dict[target][2]/exp_dict[target][0]) for target in target_genes])

    #for median_normalized
    signal_line_median = [brd4_in_2_median,brd4_in_24_median,brd4_out_2,brd4_out_24,exp_in_2,exp_in_24,exp_out_2,exp_out_24]
    #rounding
    signal_line_median = [round(x,4) for x in signal_line_median]
    new_line_median = [gene]+signal_line_median
    pca_table_median.append(new_line_median)




# out_path = '%s/tables/161202_rasmc_tf_pca.txt' % (projectFolder)
# utils.unParseTable(pca_table,out_path,'\t')

out_path_median = '%s/tables/rasmc_tf_pca_median.txt' % (projectFolder)
utils.unParseTable(pca_table_median,out_path_median,'\t')
