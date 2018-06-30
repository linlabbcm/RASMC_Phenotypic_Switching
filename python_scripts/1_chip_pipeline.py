
#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2017 YOUR NAME HERE and Charles Lin lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of YOUR PROJECT HERE




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
import itertools

from scipy import stats

#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'rasmc_all'
genome ='rn6'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files

#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Rx
chip_data_file = '%sdata_tables/RASMC_CHIP_DATA_TABLE.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I, LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)

    #assumes macs has already been run and formatted
#    run_macs(chip_data_file)

#    sys.exit()

    print('\n\n')
    print('#======================================================================')
    print('#======================II. PLOTTING BRD4 TRACKS========================')
    print('#======================================================================')
    print('\n\n')

    # for BRD4
    dataFile = chip_data_file

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_1.gff'

    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    
    
    
    plotName = 'rasmc_all_figure_1_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_2.gff'
    
    plotName = 'rasmc_all_figure_2_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')
    
    print('AAA')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3.gff'

    plotName = 'rasmc_all_figure_3_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')

    figureGFF=[['chr16','Btd','',7671879,7881204,'','.','Btd']]
    figureGFFPath='%sBtd.gff' % (gffFolder)
    utils.unParseTable(figureGFF,figureGFFPath,'\t')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Btd.gff'

    plotName = 'rasmc_all_btd_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/enrich_tfs.gff'

    plotName = 'rasmc_all_enrich_tfs_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')
    
    print('BBB')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Edn1.gff'

    plotName = 'rasmc_all_edn1_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')



    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Edn1.gff'

    plotName = 'rasmc_all_edn1_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/regions_fosl1_atf4_spry2_thbs1.gff'

    
    plotName = 'spry2_ext_brd4_tracks_with_up_down_beds'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/regions_fosl1_atf4_spry2_thbs1.gff'

    plotName = 'spry2_ext_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString,plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')


    print('CCC')



    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Acta2.gff'


    plotName = 'Acta2_brd4_tracks'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Acta2.gff'

    plotName = 'Acta2_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='',plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/cluster_b_genes_of_interest.gff'


    plotName = 'cluster_b_goi_brd4_tracks_with_up_down_beds'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/cluster_b_genes_of_interest.gff'

    plotName = 'cluster_b_goi_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')


    print('DDD')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3_supp_tracks.gff'


    plotName = 'figure_3_supp_brd4_tracks_with_up_down_beds'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3_supp_tracks.gff'

    plotName = 'figure_3_supp_pol2_tracks'
    namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW']
    print(namesList)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = bedString,plotType= 'MULTIPLE',extension=200,multiPage = False,debug=False,nameString ='' ,rpm=True,rxGenome = '')

    figureGFFPath= '/storage/cylin/grail/projects/rasmc_all/gff/thrb_enhancer_regions.gff'
    
    plotName = 'thrb_enhancer_regions_brd4'
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    groupList = []
    for name in namesList:
        if name.count('UNSTIM') > 0:
               groupList.append('UNSTIM')
        if name.count('PDGF_2H') > 0:
               groupList.append('PDGF_2H')
        if name.count('PDGF_24H') >0:
               groupList.append('PDGF_24H')

    print(namesList)
    groupString = ','.join(groupList)
    print(groupString)

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed = '',plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString ,rpm=True,rxGenome = '')


    print('EEE')




    print('\n\n')
    print('#======================================================================')
    print('#======================III. ENHANCER PROMOTER==========================')
    print('#======================================================================')
    print('\n\n')

    dataDict = pipeline_dfci.loadDataTable(chip_data_file)
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

    allLoci = []

#    for name in namesList:

#          collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#          allLoci += collection.getLoci()

    #do this for each one in the namesList
    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_BRD4_STITCHED_-0_+0.gff','\t')







    dataFile = chip_data_file
    namesList = ['RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']
    gff = '%sRN6_RASMC_BRD4_STITCHED_-0_+0.gff' % (gffFolder)
    activityTable = '%sactiveListTable.txt' % (tableFolder)
    outputFolder = '%senhancerPromoter/BRD4/' % (projectFolder)
    desc_string = '0_STITCH_-_JQ1'

#    makeEnhProBash(dataFile,gff,activityTable,namesList,outputFolder,desc_string)


    annotTable = utils.parseTable(annotFile,'\t')
    dictFile = '/storage/cylin/bin/pipeline/crc/annotation/MotifDictionary.txt'
    motifs = utils.parseTable(dictFile,'\t')
    
    TFlist = []
    
#    for line in annotTable[1:]:
#        for motif in motifs:
#            if line[12].upper()==motif[1]:
#                NMID = line[1]
#                name = line[12].upper()
#                new_line=[NMID,name]
#                TFlist.append(new_line)
    
#    utils.unParseTable(TFlist,'/storage/cylin/grail/projects/rasmc_all/tables/TFlist.txt','\t')

    print('\n\n')
    print('#======================================================================')
    print('#======================IV. TF NETWORK ANALYSIS=========================')
    print('#======================================================================')
    print('\n\n')


    #==========================================================================
    #===================ENH TF GFF=============================================
    #==========================================================================
    print('Creating Brd4 Enhancer TF in degree table from crc rasmc_h3k27ac_0_tss_ENHANCER_TF_TABLE.txt')
    
#    EnhTableFile = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_ENHANCER_TF_TABLE.txt' % (projectFolder)

#    EnhTable = utils.parseTable(EnhTableFile,'\t')
#    genesOfInterest=[]

#    enhGFF = []

    #for each region in cluster 1 (as determined by last column) add a gff line for the +/-50kb region   
    #e.g. gffLine = ['chr19','1_merged_rasmc_h3k27ac_1_lociStitched','',25834879-50000,25916534+50000,'','.','','1_merged_rasmc_h3k27ac_1_lociStitched']
#    print(EnhTable[1])

#    for line in EnhTable:
#        if line != EnhTable[0]:
#            gene = line[4]
#            ID = line[0]
#            chrom = line[1]
#            TSS = ''
#            start = line[2]
#            stop = line[3]
#            line_gff_format=[chrom, gene,ID,(int(start)-500),(int(stop)+500),'','.',TSS,gene]
#            enhGFF.append(line_gff_format)

#    enhGFF_Path = '%srasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE.gff' % (gffFolder,)
#    enhGFFTable = utils.unParseTable(enhGFF, enhGFF_Path, '\t')

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to Enhancer TF in degree gff')
#    dataFile = chip_data_file
#    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE.gff']
#    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)


    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making signal table for Brd4 Enhancer TF in degree')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_BRD4_ENH_TF_IN_DEGREE_signal_table.txt')

    #==========================================================================
    #===================EDGE TABLE GFF=========================================
    #==========================================================================
    print('Creating Edge Table GFF crc rasmc_h3k27ac_0_tss_EDGE_TABLE.txt')



#    peakTableFile = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_EDGE_TABLE.txt' % (projectFolder)

#    clusterTable = utils.parseTable(peakTableFile,'\t')
#    genesOfInterest=[]

#    clusterGFF = []

    #for each region in cluster 1 (as determined by last column) add a gff line for the +/-50kb region
    #e.g. gffLine = ['chr19','1_merged_rasmc_h3k27ac_1_lociStitched','',25834879-50000,25916534+50000,'','.','','1_merged_rasmc_h3k27ac_1_lociStitched']
#    print(clusterTable[1])

#    for line in clusterTable:
#        if line != clusterTable[0]:
#            gene = line[0]
#            ID = line[5]
#            chrom = line[2]
#            TSS = ''
#            start = line[3]
#            stop = line[4]
#            line_gff_format=[chrom, gene,ID,(int(start)-500),(int(stop)+500),'','.',TSS,'']
#            clusterGFF.append(line_gff_format)

#    clusterGFF_Path = '%srasmc_h3k27ac_0_tss_EDGE_TABLE.gff' % (gffFolder,)
#    clusterGFFTable = utils.unParseTable(clusterGFF, clusterGFF_Path, '\t')

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to Edge Table gff')
    dataFile = chip_data_file
#    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_EDGE_TABLE.gff']
#    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)


    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making brd4 signal table for Edge Table gff')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_EDGE_TABLE.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_EDGE_TABLE_Brd4_signal_table.txt')

    #==========================================================================
    #===================EDGE TABLE GFF=========================================
    #==========================================================================
    print('Creating crc subpeak Table')



    subpeakTableFile = '%scrc/rasmc_h3k27ac_0_tss/rasmc_h3k27ac_0_tss_all_subpeak.bed' % (projectFolder)

    subpeakTable = utils.parseTable(subpeakTableFile,'\t')
   
    subpeakGFF = []

    #for each region in cluster 1 (as determined by last column) add a gff line for the +/-50kb region
    #e.g. gffLine = ['chr19','1_merged_rasmc_h3k27ac_1_lociStitched','',25834879-50000,25916534+50000,'','.','','1_merged_rasmc_h3k27ac_1_lociStitched']
#    print(subpeakTable[1])

#    for line in subpeakTable:
#        gene = ''
#        ID = line[4]
#        chrom = line[0]
#        TSS = ''
#        start = line[1]
#        stop = line[2]
#        line_gff_format=[chrom, gene,ID,(int(start)-500),(int(stop)+500),'','.',TSS,'']
#        subpeakGFF.append(line_gff_format)
    
#    print(subpeakGFF[1])
#    subpeakGFF_Path = '%srasmc_h3k27ac_0_tss_subpeak_table.gff' % (gffFolder,)
#    subpeakGFFTable = utils.unParseTable(subpeakGFF, subpeakGFF_Path, '\t')

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to crc subpeak Table gff')
    dataFile = chip_data_file
    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_subpeak_table.gff']
    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)


    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making brd4 signal table for crc subpeak Table gff')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_subpeak_table.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/rasmc_h3k27ac_0_tss_subpeak_Brd4_signal_table.txt')


    print('\n\n')
    print('#======================================================================')
    print('#======================IV. POST ANALYSIS===============================')
    print('#======================================================================')
    print('\n\n')





    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to crc subpeak Table gff')
    dataFile = chip_data_file
    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff']
    brd4NamesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW',
'RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = brd4NamesList,extension=200,rpm=False)




    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making signal table for Brd4 Enhancer TF in degree')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_WCE_UNSTIM_POOLED','RASMC_WCE_BRD4_UNSTIM_NEW','RASMC_WCE_PDGF_2H_POOLED','RASMC_WCE_BRD4_PDGF_2H_NEW','RASMC_WCE_PDGF_24H_POOLED','RASMC_WCE_PDGF_24H_NEW','RASMC_BRD4_UNSTIM_REP1','RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_PDGF_2H_REP2','RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_24H_REP2','RASMC_BRD4_PDGF_24H_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/RASMC_BRD4_across_0_stitch_h3k27ac_regions.txt')


    #==========================================================================
    #====================MAP BAMS BATCH POLII==================================
    #==========================================================================
    dataFile = chip_data_file
    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/RN6_TSS_ALL_-300_+300.gff','/storage/cylin/grail/projects/rasmc_all/gff/RN6_BODY_ALL_+300_+3000.gff']
   
#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW'],extension=200,rpm=False)

    #==========================================================================
    #====================MAKE SIGNAL TABLE=====================================
    #==========================================================================
    print('Making signal table for RNA Pol2 TSS +- 300')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/RN6_TSS_ALL_-300_+300.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW'],medianNorm=False,output ='/storage/cylin/grail/projects/rasmc_all/signalTables/RASMC_POL2_RN6_TSS_ALL_-300_+300_signal_table.txt')

    print('Making signal table for RNA Pol2 gene body -300 +3000')
    dataFile = chip_data_file
    gffFile = '/storage/cylin/grail/projects/rasmc_all/gff/RN6_BODY_ALL_+300_+3000.gff'


#    pipeline_dfci.makeSignalTable(dataFile,gffFile,mappedFolder,namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW'],medianNorm=False,output='/storage/cylin/grail/projects/rasmc_all/signalTables/RASMC_POL2_RN6_BODY_ALL_+300_+3000_signal_table.txt')





   
    #==========================================================================
    #====================MAP BAMS BATCH POLII METAS============================
    #==========================================================================
    dataFile = chip_data_file
#    gffList =  ['/storage/cylin/grail/projects/rasmc_all/gff/RN6_TSS_ALL_-300_+300.gff','/storage/cylin/grail/projects/rasmc_all/gff/RN6_BODY_ALL_+300_+3000.gff']

#    pipeline_dfci.mapBamsBatch(dataFile,gffList,mappedFolder,overWrite =False,namesList = ['RASMC_POL2_UNSTIM_NEW','RASMC_POL2_PDGF_2H_NEW','RASMC_POL2_PDGF_2H_JQ1_NEW','RASMC_POL2_PDGF_24H_NEW','RASMC_POL2_PDGF_24H_JQ1_NEW'],extension=200,rpm=False)

   

    print('\n\n')
    print('#======================================================================')
    print('#======================MAKING GEO TABLES===============================')
    print('#======================================================================')
    print('\n\n')


    geoName='rasmc_chip'
    outputFolder = '/storage/cylin/grail/projects/rasmc_all/rasmc_geo/%s_geo/' % (geoName)
    namesList=[]

#    makeGEOTable(chip_data_file,wiggleFolder,macsFolder,namesList,geoName,outputFolder)



    print('\n\n')
    print('#======================================================================')
    print('#=======================XIV. DELTA OUT BY EDGE=========================')
    print('#======================================================================')
    print('\n\n')


    #calculating delta out degree by brd4 change at edges. only take edges in the top 50%
    #at least 1 dataset
    # crc_folder = '%scrc_atac/keratinocyte_combined' % (projectFolder)
    # analysis_name = 'keratinocyte_combined'
    # tf_edge_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list)


#    crc_folder = '%scrc/rasmc_h3k27ac_0_tss_all_motifs' % (projectFolder)
#    analysis_name = 'rasmc_h3k27ac_0_tss_all_motifs'
#    brd4_24_list = ['RASMC_BRD4_PDGF_24H_NEW','RASMC_BRD4_PDGF_24H_REP2']
#    brd4_2_list = ['RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_2H_REP2']
#    brd4_unstim_list=['RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_UNSTIM_REP1']
#    tf_edge_brd4_delta(crc_folder,chip_data_file,analysis_name,brd4_24_list,brd4_2_list)

    crc_folder = '%scrc/rasmc_h3k27ac_0_tss_am_ef' % (projectFolder)
    analysis_name = 'rasmc_h3k27ac_0_tss_am_ef'
#    brd4_24_list = ['RASMC_BRD4_PDGF_24H_NEW','RASMC_BRD4_PDGF_24H_REP2']
    brd4_2_list = ['RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_2H_REP2']
    brd4_unstim_list=['RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_UNSTIM_REP1']
#    tf_edge_brd4_delta(crc_folder,chip_data_file,analysis_name,brd4_2_list,brd4_unstim_list)




    print('\n\n')
    print('#======================================================================')
    print('#=======================XV. BRD4 DELTA IN =================================')
    print('#======================================================================')
    print('\n\n')


    crc_folder = '%scrc/rasmc_h3k27ac_0_tss_all_motifs' % (projectFolder)
    analysis_name = 'rasmc_h3k27ac_0_tss_all_motifs'
    brd4_24_list = ['RASMC_BRD4_PDGF_24H_NEW','RASMC_BRD4_PDGF_24H_REP2']
    brd4_2_list = ['RASMC_BRD4_PDGF_2H_NEW','RASMC_BRD4_PDGF_2H_REP2']
    brd4_unstim_list=['RASMC_BRD4_UNSTIM_NEW','RASMC_BRD4_UNSTIM_REP1']
    brd4_group_list=[brd4_unstim_list,brd4_2_list,brd4_24_list]
#    tf_brd4_in_delta(crc_folder,chip_data_file,analysis_name,brd4_group_list)



    print('\n\n')
    print('#======================================================================')
    print('#=======================XVI. CORR HEATMAP =================================')
    print('#======================================================================')
    print('\n\n')


    analysis_name = 'rasmc_h3k27ac_0_tss_all_motifs_string_db_tfs_ef'
    motifBedDir = '%scrc/rasmc_h3k27ac_0_tss_all_motifs/motif_beds/' % (projectFolder)
#    pipeline_dfci.plotCRCCorrMaps(analysis_name,motifBedDir,tf_list_path='/storage/cylin/grail/projects/rasmc_all/tables/string_db_tf_list_ef.txt',window=50)

    analysis_name = 'rasmc_h3k27ac_0_tss_all_motifs_am_ef'
    motifBedDir = '%scrc/rasmc_h3k27ac_0_tss_am_ef/motif_beds/' % (projectFolder)
#    pipeline_dfci.plotCRCCorrMaps(analysis_name,motifBedDir,tf_list_path='/storage/cylin/grail/projects/rasmc_all/tables/string_db_tf_list_ef.txt',window=50)


####################################
# '''    dataDict=pipeline_dfci.loadDataTable(chip_data_file)
#     bash_file_name='/storage/cylin/grail/projects/rasmc_all/chip_data_stats.sh'
#     bashFile=open(bash_file_name,'w')
#     bashFile.write('#!/usr/bin/bash\n')
#     names=dataDict.keys()
#     for name in names:
#         bam_name=dataDict[name]['bam']
#         cmd='samtools flagstat %s > /storage/cylin/grail/projects/rasmc_all/chip_bam_stats.txt\n' % (bam_name)
#         bashFile.write(cmd)
# '''

    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping Brd4 bams to peaks')
    dataFile = chip_data_file
    dataDict=pipeline_dfci.loadDataTable(chip_data_file)
    names=dataDict.keys()
#    for name in names:
#        if len(dataDict[name]['enrichedMacs'])>4:
#            peak_name=dataDict[name]['enrichedMacs']
#            peak_path='%s%s' % (macsEnrichedFolder,peak_name)
#            gff_path='%s%s.gff' % (gffFolder,peak_name.split('.bed')[0])
#            utils.bedToGFF(peak_path,output=gff_path)
#            gffList=[gff_path]
#            namesL=[name]
#            pipeline_dfci.mapBamsBatch(dataFile, gffList,mappedFolder,overWrite=True,namesList=namesL,extension=0,rpm=True)


    namesL=names
    tss_gff_path='%sRN6_TSS_ALL_-300_+300.gff' % (gffFolder)
    gffList=[tss_gff_path]
    pipeline_dfci.mapBamsBatch(dataFile, gffList,mappedFolder,overWrite=False,namesList=namesL,extension=0,rpm=True)

        
       
#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RUNNING MACS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def run_macs(dataFile):
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    namesList.sort()
    print(namesList)
    pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9')
    os.chdir(projectFolder) # the silly call macs script has to change into the output dir
    #so this takes us back to the project folder

    #to check for completeness, we will try to find all of the peak files
    peak_calling_done = False
    while not peak_calling_done:
        dataDict = pipeline_dfci.loadDataTable(dataFile)
        namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
        for name in namesList:
            peak_path = '%s%s/%s_summits.bed' % (macsFolder,name,name)
            print('searching for %s' % (peak_path))
            if utils.checkOutput(peak_path,1,180):
                peak_calling_done =True
                print('found %s' % (peak_path))
                continue
            else:
                print('Error: peak calling timed out')
                sys.exit()
    
    #now format the macs output
    print('formatting macs output')
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink ='',useBackground=True)
    print('Finished running Macs 1.4.2')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAKING ENHANCER PROMOTER BASH SCRIPTS~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def makeEnhProBash(dataFile,gff='',activityTable='',namesList=[],outputFolder='',desc_string=''):
    dataDict = pipeline_dfci.loadDataTable(dataFile)

    print(dataDict)
    if len(namesList) == 0:
        namesList = dataDict.keys()
    namesList.sort()
    print(namesList)
    
    for name in namesList:
        bashFileName = "%s%s_%s_enhPro.sh" % (outputFolder,name,desc_string)
        print('Writing %s') % bashFileName
        bashFile = open(bashFileName,'w')

        #shebang
        bashFile.write('#!/usr/bin/bash\n')

        line = 'BAMS=\'%s\'' % (dataDict[name]['bam'])
        bashFile.write(line+'\n')

        line = 'CONTROLS=\'%s\'' % (dataDict[name]['background'])
        bashFile.write(line+'\n\n')

        line = 'OUTPUT=\'./\''
        bashFile.write(line+'\n')

        line = 'INPUT=\'%s\'' % (gff)
        bashFile.write(line+'\n')
        line = 'ACTIVITY=\'%s\'' % (activityTable)
        bashFile.write(line+'\n')
        line = 'NAME=\'%s_%s\'' % (name,desc_string)
        bashFile.write(line+'\n\n\n')

        line = 'python /storage/cylin/bin/pipeline/enhancerPromoter.py -b $BAMS -g %s -i $INPUT -o $OUTPUT -a $ACTIVITY --name $NAME' % (dataDict[name]['genome'])
        bashFile.write(line+'\n')
        
        bashFile.close()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~MAKING GEO TABLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeGEOTable(dataFile,wiggleFolder,macsFolder,namesList,geoName,outputFolder =''):

    '''
    makes a geo table and a bash script to format everything
    '''
    dataDict = pipeline_dfci.loadDataTable(dataFile)

    #first make a reverse wce dict

    backgroundDict = {}
    if len(namesList) == 0:
        namesList = dataDict.keys()
    for name in namesList:

        background = dataDict[name]['background']
        backgroundDict[background] = name


    outputFolder = pipeline_dfci.formatFolder(outputFolder,True)
    bashFileName = '%s%s_bash.sh' % (outputFolder,geoName)

    bashFile = open(bashFileName,'w')

    geoTable = [['SAMPLE_NAME','TITLE','CELL_TYPE','PROCESSED_FILE','RAW_FILE','BARCODE']]

    namesList.sort()

    for name in namesList:

        print(name)

        sampleName = dataDict[name]['uniqueID']
        title = name
        cell_type = name.split('_')[0]
        processed_file = "%s.wig.gz" % (name)
        raw_file = "%s.fastq.gz" % (name)
        
        fastqFile = dataDict[name]['fastq']
        uniqueID = dataDict[name]['uniqueID']
        try:
            barcode = pipeline_dfci.getTONYInfo(uniqueID,38)
        except IndexError:
            barcode = ''
        newLine = [sampleName,title,cell_type,processed_file,raw_file,barcode]
        print(newLine)
        geoTable.append(newLine)


    utils.unParseTable(geoTable,"%s%s_meta.xls" % (outputFolder,geoName),'\t')

    #now make the folder to hold everything and the relevant bash script
    if len(outputFolder) == 0:
        outputFolder ='./%s/' % (geoName)

    else:
        outputFolder = outputFolder + geoName + '/'

    pipeline_dfci.formatFolder(outputFolder,True)

    wiggleFolder = pipeline_dfci.formatFolder(wiggleFolder,False)
    macsFolder = pipeline_dfci.formatFolder(macsFolder,False)

    #now make the bash file
    bashFile.write('#!/usr/bin/bash\n')
    bashFile.write('cd %s\n' %(outputFolder))
    bashFile.write('\n')

    #write the untar command
    for name in namesList:
        fastqFile = dataDict[name]['fastq']
        if len(fastqFile) == 0:
            print "NO FASTQ FILE FOR %s" % (name)
            continue
        if fastqFile.count('tar.gz') > 0: #for files generated by whitehead that have tar header #####RACHEL READ HERE
            tarCmd = 'tar --strip-components 5 --to-stdout -xzvf %s | gzip > %s.fastq.gz\n' % (fastqFile,name)
        else:
            tarCmd = 'cp %s %s.fastq.gz\n' % (fastqFile,name)
        bashFile.write(tarCmd)

    bashFile.write('\n\n\n')
    #write the wiggle cp command
    for name in namesList:
        if name.count('WCE') == 1 or name.count('INPUT') == 1 :
            refName = backgroundDict[name]
            controlWiggleFile = '%s%s/%s_MACS_wiggle/control/%s_control_afterfiting_all.wig.gz' % (macsFolder,refName,refName,refName)
            wigCmd = "#cp '%s' %s.wig.gz\n" % (controlWiggleFile,name)
            #wigCmd = "cp '%swceWiggles/%s_control_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,refName,name)
        else:
            wigCmd = "#cp '%s%s_treat_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,name,name)
        bashFile.write(wigCmd)

    #write the md5sums for the wiggles
    bashFile.write('\n\n\n')
    bashFile.write("echo '' > md5sum.txt\n")
    for name in namesList:
        md5Cmd = 'md5sum %s.wig.gz >> md5sum.txt\n' % (name)
        bashFile.write(md5Cmd)

    #write md5sums for the fastqs
    for name in namesList:
        md5Cmd = 'md5sum %s.fastq.gz >> md5sum.txt\n' % (name)
        bashFile.write(md5Cmd)

    #the big tar command
    tarCmd = '#tar -cvzf %s.tar.gz %s\n' % (geoName,outputFolder)
    bashFile.write(tarCmd)
    bashFile.close()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF EDGES~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list,output=''):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence this is by subpeaks
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    edge_path = '%s%s_EDGE_TABLE.txt' % (crc_folder,analysis_name)

    #make a gff of the edge table
    edge_table = utils.parseTable(edge_path,'\t')
    edge_gff = []
    for line in edge_table[1:]:
        gff_line = [line[2],'%s_%s' % (line[0],line[1]),'',line[3],line[4],'','.','','%s_%s' % (line[0],line[1])]
        edge_gff.append(gff_line)

    edge_gff_path = '%s%s_EDGE_TABLE.gff' % (crc_folder,analysis_name)
    utils.unParseTable(edge_gff,edge_gff_path,'\t')


    #direct the output to the crc folder
    signal_path = '%s%s_EDGE_TABLE_signal.txt' % (crc_folder,analysis_name)



    all_brd4_list = y_brd4_list + o_brd4_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[edge_gff_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))



    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 vs 0  signal table at edges')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)

    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]
    o_signal_vector = []
    y_signal_vector = []
    for line in signal_table[1:]:
        o_signal = numpy.mean([float(line[col]) for col in o_columns])
        y_signal = numpy.mean([float(line[col]) for col in y_columns])

        o_signal_vector.append(o_signal)
        y_signal_vector.append(y_signal)

    o_median = numpy.median(o_signal_vector)
    y_median = numpy.median(y_signal_vector)

    print('old median BRD4 signal')
    print(o_median)
    print('young median BRD4 signal')
    print(y_median)

    #now that we have the median, we can take edges where at least 1 edge is above the median
    #and both are above zero and generate a new table w/ the fold change

    signal_filtered_path = string.replace(signal_path,'.txt','_filtered.txt')
    if utils.checkOutput(signal_filtered_path,0,0):
        print('Found filtered signal table for edges at %s' % (signal_filtered_path))
        signal_table_filtered = utils.parseTable(signal_filtered_path,'\t')
    else:

        signal_table_filtered = [signal_table[0]+['O_BRD4_MEAN','Y_BRD4_MEAN','Y_vs_O_LOG2']]
        for line in signal_table[1:]:
            o_signal = numpy.mean([float(line[col]) for col in o_columns])
            y_signal = numpy.mean([float(line[col]) for col in y_columns])

            if (o_signal > o_median or y_signal > y_median) and min(o_signal,y_signal) >0:
                delta = numpy.log2(y_signal/o_signal)
                new_line = line + [o_signal,y_signal,delta]
                signal_table_filtered.append(new_line)

        utils.unParseTable(signal_table_filtered,signal_filtered_path,'\t')

    #now get a list of all TFs in the system
    tf_list = utils.uniquify([line[0].split('_')[0] for line in signal_table_filtered[1:]])
    tf_list.sort()
    print(tf_list)

    out_degree_table = [['TF_NAME','EDGE_COUNT','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD','DELTA_SEM']]

    for tf_name in tf_list:
        print(tf_name)
        edge_vector = [float(line[-1]) for line in signal_table_filtered[1:] if line[0].split('_')[0] == tf_name]

        edge_count = len(edge_vector)
        delta_mean = round(numpy.mean(edge_vector),4)
        delta_median = round(numpy.median(edge_vector),4)
        delta_std = round(numpy.std(edge_vector),4)
        delta_sem = round(stats.sem(edge_vector),4)
        tf_out_line = [tf_name,edge_count,delta_mean,delta_median,delta_std,delta_sem]
        out_degree_table.append(tf_out_line)

    if output == '':
        #set final output
        output_path = '%s%s_EDGE_DELTA_OUT.txt' % (crc_folder,analysis_name)

    else:
        output_path = output

    utils.unParseTable(out_degree_table,output_path,'\t')
    print(output_path)

    return(output_path)



#=========================================================================
#====================DELTA BRD4 IN TABLE==================================
#=========================================================================

def tf_brd4_in_delta(crc_folder,chip_dataFile,analysis_name,brd4_group_list):
    '''
    calculates changes in brd4 in degree for tfs from crc output
    brd_group_list should be ranked from smallest timepoint to greatest
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    enh_tf_path = '{}{}{}'.format(crc_folder,analysis_name,'_ENHANCER_TF_TABLE.txt')

    #make gff of enhancer tf table
    enh_tf_table = utils.parseTable(enh_tf_path,'\t')
    enh_gff = []
    for line in enh_tf_table[1:]:
        gff_line = [line[1],line[4],line[0],(int(line[2])-500),(int(line[3])+500),'','.','',line[4]]
        enh_gff.append(gff_line)

    enh_tf_gff_path = '{}{}{}'.format(crc_folder,analysis_name,'_ENH_TF_TABLE.gff')
    utils.unParseTable(enh_gff,enh_tf_gff_path,'\t')

    #direct output to the crc folder
    signal_path = '{}{}{}'.format(crc_folder,analysis_name,'_ENH_TF_TABLE_signal.txt')

    all_brd4_list = []
    for group in brd4_group_list:
        for item in group:
            all_brd4_list.append(item)

    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[enh_tf_gff_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('{}{}'.format('Found previous signal table at ',signal_path))


    #now bring in the signal table as a dictionary using the tf name as the id
    print('making log2 fc signal table of enhnacer tfs')
    signal_table = utils.parseTable(signal_path, '\t')
    signal_dict = defaultdict(float)

    medians = []

    for group in brd4_group_list:
        group_cols = [signal_table[0].index(name) for name in group]
        group_signal_vector =[ ]
        for line in signal_table[1:]:
            group_signal = numpy.mean([float(line[col]) for col in group_cols])

            group_signal_vector.append(group_signal)
        med = numpy.median(group_signal_vector)
        medians.append(med)

    print(medians)

    delta_in_path = '{}{}{}'.format(crc_folder,analysis_name,'_TF_DELTA_IN.txt')

    group_names = []
    for group in brd4_group_list:
        name_split = group[0].split('_')
        name = '_'.join(name_split[:-1])
        print(name)
        group_names.append(name)
    print(group_names)
    combos = list(itertools.combinations(group_names,2))
    log_names = []
    for combo in combos:
        log_name = '{}{}{}{}'.format(combo[1],'_vs_',combo[0],'_LOG2')
        log_names.append(log_name)
    foo='GENE_ID'
    delta_in_table=[]
    delta_in_table.append(foo)
    delta_in_table = [delta_in_table+group_names+log_names]
    print(delta_in_table)
    groups = len(brd4_group_list)
#    print(groups)
    iter_list = list(range(1,groups+1))
#    print(iter_list)
    col_combos = list(itertools.combinations(iter_list,2))
#    print(col_combos)
    count=0
    for line in signal_table:
#        print(count)
        count=count+1
        if line == signal_table[0]:
            header=signal_table[0]
        if line != signal_table[0]:
            new_line=[line[0]]
            for group in brd4_group_list:
                group_cols = [header.index(name) for name in group]
                group_signal = numpy.mean([float(line[col]) for col in group_cols])
                new_line.append(group_signal)
            log_line=[]
            #print(new_line)
            for combo in col_combos:
                a=new_line[combo[1]]
                b=new_line[combo[0]]
                #print(a)
                #print(b)
                log_fc = numpy.log2(a/b)
                #print(log_fc)
                log_line.append(log_fc)
                #print(log_line)
            new_line = new_line + log_line
            #print(new_line)
            delta_in_table.append(new_line)

    del_table=[delta_in_table[0]]
    for line in delta_in_table[1:]:
        name_line=[]
        name=line[0]
        print(name)
        length = len(line)
        split_name = name.split(',')
        if len(split_name) > 1:
        #    print(split_name)
            for i in split_name:
                name_line.append(i)
                name_line=name_line+line[1:(length-1)]
                print(name_line)
                del_table.append(name_line)
                name_line=[]
        else:
            name_line = line
            del_table.append(name_line)


    utils.unParseTable(del_table,delta_in_path,'\t')
    print(delta_in_path)
    return(delta_in_path)



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
