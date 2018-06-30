
import sys
sys.path.append('/storage/cylin/bin/pipeline/')

import pipeline_dfci
import utils
import string

from collections import defaultdict

#Paths to tables and files
annotFile = '/storage/cylin/bin/pipeline/annotation/rn6_refseq.ucsc'
enhTablePath = '/storage/cylin/grail/projects/rasmc_all/crc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_ENHANCER_TABLE.txt'
expFile = '/storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_exprs_fpkm_means.txt'
edgeFile = '/storage/cylin/grail/projects/rasmc_all/TF_Network/filtered_edge_table.txt'

#Parse Tables
enhTable = utils.parseTable(enhTablePath, '\t')
expTable = utils.parseTable(expFile, '\t')
edgeTable = utils.parseTable(edgeFile,'\t')

geneList = []

print(enhTable[1])

#Create startDict for enhancer table genes
startDict = utils.makeStartDict(annotFile)

#print(startDict.keys())

#Create expDict to enforce >10 FKPM cutoff in one of 0H,2H+PDGF,24H+PDGF
expDict = defaultdict(list)

for line in expTable[1:]:
    expDict[line[0].upper()] = [float(line[1]),float(line[2]),float(line[4])]

print(len(expDict))

def getClosest(enhancer_table, ref_dictionary):
    closestGenes = []
    j = 0
    for line in enhancer_table[1:10]:
        print(j)
        j = j+1
        regionID = line[0]
        #print(regionID)
        tss_enh = int(line[2])
        #print(tss_enh)
        if len(line) > 4:
            target_genes = line[4].split(',')
        else: 
            closest_gene = ''
   
        if len(target_genes) > 0:
            refIDList = utils.nameToRefseq(target_genes,annotFile,unique=False)
            for i in range(len(refIDList)):
                newRefID = refIDList[i]
                newTSS = startDict[newRefID[1]]['start']
               # print(newRefID)
                if i == 0:
                    tss_index = 0
                else:
                    tss_diff = abs(int(newTSS[0]) - tss_enh)
                    old_diff = abs(int(oldTSS[0]) - tss_enh)
                    if tss_diff < old_diff:
                        tss_index = i
                
                refID = refIDList[tss_index]
                oldTSS = startDict[refID[1]]['start']
                #print(refID)
                closest_gene = startDict[refID[1]]['name'].upper() 
                print(closest_gene)
                
        gene_line = [regionID,closest_gene]
        closestGenes.append(gene_line)
        
    closestTable = utils.unParseTable(closestGenes, '/storage/cylin/grail/projects/rasmc_all/TF_Network/enh_table_closest_genes.txt', '\t')
    return closestGenes


#getClosest(enhTable, startDict)

#print(closestGenes[0:9])
 
closest_genes = utils.parseTable('/storage/cylin/grail/projects/rasmc_all/TF_Network/enh_table_closest_genes.txt','\t')

print(expDict[1])

enhDict = defaultdict(list)

for line in closest_genes:
    if len(line)>1:
        gene = line[1]
    else: 
        gene = ''
    enhDict[line[0]]=[gene]

filtered_edge = []

for edge in edgeTable:
    print(edge[5])
    if len(edge[5]) > 0:
        closest_gene = (enhDict[edge[5]][0]).upper()
        print(closest_gene)
        if len(closest_gene) > 0:
            if expDict[closest_gene][0] >= 10 or expDict[closest_gene][1] >= 10 or expDict[closest_gene][2] >= 10:
                new_line = [edge[0],closest_gene,edge[2],edge[3],edge[4],edge[5],edge[6]]
                #print(new_line)
                filtered_edge.append(new_line)

utils.unParseTable(filtered_edge,'/storage/cylin/grail/projects/rasmc_all/TF_Network/filtered_edge_exp_closest_genes_filter.txt','\t')



