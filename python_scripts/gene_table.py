import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils

import argparse
import cPickle
import utils
import pipeline_dfci
import subprocess
import os
import string
import tempfile
import zlib
import numpy
import re
import time
from distutils.spawn import find_executable
from collections import defaultdict


def makeGeneTable(peakTable,analysisName):

    '''
    takes the peak table and makes a gene centric table

    '''


    geneDict = {}

    geneTable = [['GENE','%s_TSS_SIGNAL' % (analysisName),'%s_DISTAL_SIGNAL' % (analysisName),'regionID']]

    #now iterate through the table
    for line in peakTable[1:]:
        regionLength = int(line[4])

        signal = float(line[9]) * regionLength


        #genes where this particular peak overlaps the tss 1kb window
        #where there are both overlap and proximal meet
        if len(line) == 15:
            overlapGeneList = [gene for gene in line[-2].split(',') if len(gene) > 0]
            if overlapGeneList.count('107'):
                print(line)
                sys.exit()
            for overlapGene in overlapGeneList:
                if geneDict.has_key(overlapGene) == False:
                    geneDict[overlapGene] = {'tss':0.0,'distal':0.0,'regionID':[]}
                    geneDict[overlapGene]['regionID'] += [line[0]]
                #there can be a nasty 1 overlap case where the region might overlap by the overlapping gene list, but not be real
                if int(line[5]) == 1:
                    geneDict[overlapGene]['tss'] += signal
                    geneDict[overlapGene]['regionID'] += [line[0]]
                else: #this is the case where the mycn site is just outside of the promoter or overlapping the gene locus/body these are rar
                    geneDict[overlapGene]['distal'] += signal
                    geneDict[overlapGene]['regionID'] += [line[0]]

            proximalGeneList = [gene for gene in line[-1].split(',') if len(gene) > 0]
            for proximalGene in proximalGeneList:
                if geneDict.has_key(proximalGene) == False:
                    geneDict[proximalGene] = {'tss':0.0,'distal':0.0,'regionID':[]}
                    geneDict[proximalGene]['regionID'] += [line[0]]
                if int(line[5]) == 0:
                    geneDict[proximalGene]['distal'] += signal
                    geneDict[proximalGene]['regionID'] += [line[0]]
        #where there's just overlap
        if len(line) == 14:
            overlapGeneList = [gene for gene in line[-1].split(',') if len(gene) > 0]
            if overlapGeneList.count('107'):
                print(line)
                sys.exit()
            for overlapGene in overlapGeneList:
                if geneDict.has_key(overlapGene) == False:
                    geneDict[overlapGene] = {'tss':0.0,'distal':0.0,'regionID':[]}
                    geneDict[overlapGene]['regionID'] += [line[0]]
                #there can be a nasty 1 overlap case where the region might overlap by the overlapping gene list, but not be real
                if int(line[5]) == 1:
                    geneDict[overlapGene]['tss'] += signal
                    geneDict[overlapGene]['regionID'] += [line[0]]
                else: #this is the case where the mycn site is just outside of the promoter or overlapping the gene locus/body these are rar
                    geneDict[overlapGene]['distal'] += signal
                    geneDict[overlapGene]['regionID'] += [line[0]]



    geneList = geneDict.keys()
    geneList.sort()

    for gene in geneList:
        newLine = [gene]
        newLine.append(geneDict[gene]['tss'])
        newLine.append(geneDict[gene]['distal'])
        newLine.append(','.join(geneDict[gene]['regionID']))
        geneTable.append(newLine)

    utils.unParseTable(geneTable,'/storage/cylin/grail/projects/rasmc_all/tables/gene_to_regions.txt','\t')
    return geneTable


peak_path = '/storage/cylin/grail/projects/rasmc_all/enhancerPromoter/BRD4_H3K_GFF/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_h3k_gff/RASMC_BRD4_PDGF_2H_NEW_0_STITCH_-_JQ1_h3k_gff_PEAK_TABLE.txt'
peakTable = utils.parseTable(peak_path,'\t')

makeGeneTable(peakTable,'test')


