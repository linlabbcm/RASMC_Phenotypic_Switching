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
chiprx_data_file = '%sdata_tables/RASMC_H3K27AC_CHIPRX_ALL.txt' % (projectFolder)




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
    pipeline_dfci.summary(chiprx_data_file)

    #assumes macs has already been run and formatted
    #run_macs(chiprx_data_file)

    #sys.exit()
    print('\n\n')
    print('#======================================================================')
    print('#=====================II. MAPPING H3K27AC REGIONS======================')
    print('#======================================================================')
    print('\n\n')

    #assumes macs has already been run and formatted
    #run_macs(chiprx_data_file)
    
    #next makes gffs of the k27ac regions individually, plus of intersect and union
    set_name = 'RASMC_H3K27AC'
    names_list = ['RASMC_H3K27AC_UNSTIM_NEW',
                  'RASMC_H3K27AC_PDGF_2H_NEW',
                  'RASMC_H3K27AC_PDGF_24H_NEW',
                  'RASMC_H3K27AC_UNSTIM_REP1',
                  'RASMC_H3K27AC_PDGF_2H_REP2',
                  'RASMC_H3K27AC_PDGF_24H_REP2',
                  ]

    #gff_list= makeStitchedGFF(chiprx_data_file,set_name,names_list,write_individual=True)
    #print(gff_list)
    #sys.exit()

    gff_list =  ['%sRN6_RASMC_H3K27AC_UNSTIM_NEW_-0_+0.gff' % (gffFolder),
                 '%sRN6_RASMC_H3K27AC_UNSTIM_REP1_-0_+0.gff' % (gffFolder),
                 '%sRN6_RASMC_H3K27AC_INTERSECT_-0_+0.gff' % (gffFolder),
                 '%sRN6_RASMC_H3K27AC_UNION_-0_+0.gff' % (gffFolder),
                ]
    #print(gff_list)
    #signal_table_list=pipeline_dfci.map_regions(chiprx_data_file,gff_list,mappedFolder,signalFolder,names_list=[],medianNorm=False,output='')
    #print(signal_table_list)
    
    #sys.exit()
    
    print('\n\n')
    print('#======================================================================')
    print('#================III. CALCULATING CHIPRX SCALE FACTORS=================')
    print('#======================================================================')
    print('\n\n')

    scale_table_path = '%sRN6_RASMC_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)
    #scale_table_table = writeScaleFactors(chiprx_data_file,namesList=[],primary_genome='rn6',rx_genome='dm6',output=scale_table_path)

    #sys.exit()

    print('\n\n')
    print('#======================================================================')
    print('#======================IV. MAKING BOX PLOTS============================')
    print('#======================================================================')
    print('\n\n')



    #making boxplots
    boxplot_script_path = '%sr_scripts/4_chiprx_plots.R' % (projectFolder)
    output_folder = utils.formatFolder('%sfigures/4_chiprx_plots/' % (projectFolder),True)
    scale_table_path = '%sRN6_RASMC_NEW_CHIPRX_SCALE_FACTORS.txt' % (tableFolder)
    #=============================================================================

    #for new unstim h3k27ac regions    
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_UNSTIM_NEW'
    signal_table_path = '%sRN6_RASMC_H3K27AC_UNSTIM_NEW_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)
    
    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)

    #for old unstim h3k27ac regions
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_UNSTIM'
    signal_table_path = '%sRN6_RASMC_H3K27AC_UNSTIM_REP1_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)

    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)


    #for rasmc intersect regions
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_INTERSECT'
    signal_table_path = '%sRN6_RASMC_H3K27AC_INTERSECT_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)
    
    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)


    #for rasmc union regions
    names_string = 'RASMC_H3K27AC_UNSTIM_REP1,RASMC_H3K27AC_UNSTIM_NEW,RASMC_H3K27AC_PDGF_2H_REP2,RASMC_H3K27AC_PDGF_2H_NEW,RASMC_H3K27AC_PDGF_24H_REP2,RASMC_H3K27AC_PDGF_24H_NEW'
    background_string = 'RASMC_WCE_UNSTIM_REP2,RASMC_WCE_UNSTIM_NEW,RASMC_WCE_PDGF_2H_REP2,RASMC_WCE_PDGF_2H_NEW,RASMC_WCE_PDGF_24H_REP2,RASMC_WCE_PDGF_24H_NEW'
    rep_string = 'RASMC_H3K27AC_UNSTIM,RASMC_H3K27AC_PDGF_2H,RASMC_H3K27AC_PDGF_24H'
    plot_name = 'RASMC_H3K27AC_UNION'
    signal_table_path = '%sRN6_RASMC_H3K27AC_UNION_-0_+0_RASMC_H3K27AC_CHIPRX_ALL_SIGNAL.txt' % (signalFolder)
    
    r_cmd = 'Rscript %s %s %s %s %s %s %s %s' % (boxplot_script_path,signal_table_path,scale_table_path,names_string,background_string,plot_name,projectFolder,rep_string)
    #os.system(r_cmd)
    #print(r_cmd)

#   sys.exit()


    print('\n\n')
    print('#======================================================================')
    print('#======================V. MAKING RN6 GFFS=============================')
    print('#======================================================================')
    print('\n\n')


#    pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species='RN6')


    print('\n\n')
    print('#======================================================================')
    print('#======================VI. MAPPING -JQ1 BAMS TO TSS REGIONS============')
    print('#======================================================================')
    print('\n\n')
    
    dataFile = chiprx_data_file
    gffList = ['/storage/cylin/grail/projects/rasmc_all/gff/RN6_TSS_ALL_-1000_+1000.gff']
    setName = 'RASMC_H3K27AC_-JQ1'
    cellTypeList = ['RASMC']
    enrichedFolder = macsEnrichedFolder
#    namesList=['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW','RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW']

#    pipeline_dfci.mapEnrichedToGFF(dataFile,setName,gffList,cellTypeList,enrichedFolder,mappedEnrichedFolder,macs=True,namesList=['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW','RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW'],useBackground=True)


#    setList = [['RASMC_H3K27AC_UNSTIM_REP1'],['RASMC_H3K27AC_UNSTIM_NEW'],['RASMC_H3K27AC_PDGF_2H_REP2'],['RASMC_H3K27AC_PDGF_2H_NEW'], ['RASMC_H3K27AC_PDGF_24H_REP2'],['RASMC_H3K27AC_PDGF_24H_NEW']]
#    mappedEnrichedFile = '%smappedEnriched/RN6_TSS_ALL_-1000_+1000/RN6_TSS_ALL_-1000_+1000_RASMC_H3K27AC_-JQ1.txt' % (projectFolder)
#    output = '%stables/RN6_H3K27AC_-JQ1_TSS_ACTIVE_-1000_+1000.txt' % (projectFolder)
#    pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)






    print('\n\n')
    print('#======================================================================')
    print('#======================VII. MAKING ACTIVITY TABLE=======================')
    print('#======================================================================')
    print('\n\n')
    
    h3k27ac_mapped_path = '%stables/RN6_H3K27AC_-JQ1_TSS_ACTIVE_-1000_+1000.txt' % (projectFolder)
    expr_path = '%scufflinks/rasmc_rna_cuffnorm/output/rasmc_rna_exprs_fpkm_means.txt' % (projectFolder)
    
    activeTable = makeActiveTable(h3k27ac_mapped_path,expr_path)


    print('\n\n')
    print('#======================================================================')
    print('#======================VIII. MAKING H3K27AC TRACKS=====================')
    print('#======================================================================')
    print('\n\n')





#    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_1.gff'
#    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_2.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


#    plotName = 'rasmc_all_figure_1_h3k27ac_tracks_merge'
    plotName= 'rasmc_all_figure_2_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'
    
    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'
    
#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString =scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


    plotName= 'rasmc_all_figure_3_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)


    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Btd.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


    plotName= 'rasmc_all_btd_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)




    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Edn1.gff'

    # # # #for H3K27AC
    dataFile = chiprx_data_file


    plotName= 'rasmc_all_edn1_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/Acta2.gff'

    dataFile = chiprx_data_file


    plotName= 'rasmc_all_Acta2_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed ='', plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/spry2_20kb_ext.gff'

    dataFile = chiprx_data_file


    plotName= 'spry2_ext_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'


#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)




    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/cluster_b_genes_of_interest.gff'

    dataFile = chiprx_data_file


    plotName= 'cluster_b_goi_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'


#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/figure_3_supp_tracks.gff'

    dataFile = chiprx_data_file


    plotName= 'figure_3_supp_h3k27ac_tracks_merge'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]
    groupString = 'UNSTIM,UNSTIM,PDGF_2H,PDGF_2H,PDGF_JQ1_2H,PDGF_JQ1_2H,PDGF_24H,PDGF_24H,PDGF_JQ1_24H,PDGF_JQ1_24H'

    scale_string = '0.5735,2.867,2.1121,1.9088,1.7726,2.4443,2.7222,1.5097,1.2819,1.9499'

    bedString = '/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_up.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v2_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_0v24_down.bed,/storage/cylin/grail/projects/rasmc_all/beds/enhPro_h3k_gff_regions_BRD4_2v24_down.bed'


#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)

    figureGFFPath = '/storage/cylin/grail/projects/rasmc_all/gff/stitching_evidence_Rara.gff'

    dataFile = chiprx_data_file


    plotName= 'enhancer_stitching_comp_peaks'
    outputFolder = utils.formatFolder('%sgenePlot' % (projectFolder),True)
    namesList = ['RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',]


    groupString = 'PDGF_24H,PDGF_24H'

    scale_string = '2.7222,1.5097'

    bedString='/storage/cylin/grail/projects/rasmc_all/macsEnriched/RASMC_H3K27AC_PDGF_24H_REP2_peaks.bed,/storage/cylin/grail/projects/rasmc_all/rose/RASMC_H3K27AC_PDGF_24H_REP2_ROSE/RASMC_H3K27AC_PDGF_24H_REP2_peaks_Enhancers_withSuper.bed'


#    pipeline_dfci.callBatchPlot(dataFile,figureGFFPath,plotName,outputFolder,namesList,uniform=True,bed =bedString, plotType= 'MERGE',extension=200,multiPage = False,debug=False,nameString = groupString,rpm=False,rxGenome='',scaleFactorString=scale_string)



    print('\n\n')
    print('#======================================================================')
    print('#======================IX. MAKING 0 STITCHED GFFS======================')
    print('#======================================================================')
    print('\n\n')


#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW','RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

    #do this for each one in the namesList
    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff','\t')

#######################Making Timepoint GFFs############################################

    ##Unstim##

#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_UNSTIM_REP1','RASMC_H3K27AC_UNSTIM_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

#    #do this for each one in the namesList
#    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_UNSTIM_NO_JQ1_STITCHED_-0_+0.gff','\t')

    ##2H+PDGF##

#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_PDGF_2H_REP2','RASMC_H3K27AC_PDGF_2H_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

#    #do this for each one in the namesList
#    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_2H_NO_JQ1_STITCHED_-0_+0.gff','\t')


    ##24H+PDGF##

#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)
#    namesList = ['RASMC_H3K27AC_PDGF_24H_REP2','RASMC_H3K27AC_PDGF_24H_NEW']

#    allLoci = []
#    for name in namesList:

#        collection = utils.importBoundRegion('/storage/cylin/grail/projects/rasmc_all/macsEnriched/%s_peaks.bed' %(name),name)

#        allLoci += collection.getLoci()

#    #do this for each one in the namesList
#    #then make a giant collection

#    giant_collection = utils.LocusCollection(allLoci,50)

#    stitched_collection = giant_collection.stitchCollection()

#    gff = utils.locusCollectionToGFF(stitched_collection)

#    utils.unParseTable(gff,'/storage/cylin/grail/projects/rasmc_all/gff/RN6_RASMC_H3K27AC_24H_NO_JQ1_STITCHED_-0_+0.gff','\t')



    print('\n\n')
    print('#======================================================================')
    print('#======================X. ENHANCER PROMOTER============================')
    print('#======================================================================')
    print('\n\n')

    dataFile = chiprx_data_file
    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]

    gff = '%sRN6_RASMC_H3K27AC_NO_JQ1_STITCHED_-0_+0.gff' % (gffFolder)
    activityTable = '%sactiveListTable.txt' % (tableFolder)
    outputFolder = '%senhancerPromoter/H3K27AC/' % (projectFolder)
    desc_string = '0_STITCH_-_JQ1'

#    makeEnhProBash(dataFile,gff,activityTable,namesList,outputFolder,desc_string)


    print('\n\n')
    print('#======================================================================')
    print('#======================XI. CALLING ROSE2===============================')
    print('#======================================================================')
    print('\n\n')

    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
                 'RASMC_H3K27AC_UNSTIM_NEW',
                 'RASMC_H3K27AC_PDGF_2H_REP2',
                 'RASMC_H3K27AC_PDGF_2H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_2H_NEW',
                 'RASMC_H3K27AC_PDGF_24H_REP2',
                 'RASMC_H3K27AC_PDGF_24H_NEW',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_REP2',
                 'RASMC_H3K27AC_PDGF_JQ1_24H_NEW'
                 ]

#    parentFolder = utils.formatFolder('%srose' % (projectFolder),True)
#    bashFileName = '%sRASMC_ENHANCER_H3K27AC_0_tss_rose.sh' % (parentFolder)


#    pipeline_dfci.callRose2(chiprx_data_file,macsEnrichedFolder,parentFolder,namesList,extraMap = [],inputFile='',tss=0,stitch='',bashFileName =bashFileName,mask='')


    print('\n\n')
    print('#======================================================================')
    print('#======================XII. CLUSTERING ROSE OUTPUT======================')
    print('#======================================================================')
    print('\n\n')




    #load the datadict
#    dataDict = pipeline_dfci.loadDataTable(chiprx_data_file)

    #analysis name
#    analysisName = 'rasmc_h3k27ac_0_tss'

    #nameslist
#    namesList = ['RASMC_H3K27AC_UNSTIM_REP1',
#                 'RASMC_H3K27AC_UNSTIM_NEW',
#                 'RASMC_H3K27AC_PDGF_2H_REP2',
#                 'RASMC_H3K27AC_PDGF_2H_NEW',
#                 'RASMC_H3K27AC_PDGF_24H_REP2',
#                 'RASMC_H3K27AC_PDGF_24H_NEW',
#                 ]

#    namesString = string.join(namesList,',')


#    print namesList


    #set up the output
#    outputFolder = '%sclustering/%s/' % (projectFolder,analysisName)
#    outputFolder = utils.formatFolder(outputFolder,True)

    #get the rose folder
#    roseFolder ='%srose/' % (projectFolder)

    #set up the bash file
#    bashFileName = '%s%s_0_tss_clustering.sh' % (outputFolder,analysisName)
#    bashFile = open(bashFileName,'w')

#    bashFile.write('#!/usr/bin/bash\n')

    #for now change into pipelinedir just to be safe
#    bashFile.write('cd /storage/cylin/bin/pipeline/\n')

    #if you want to do w/ all do this
    #this will run w/ default parameters
    #see the documentation for available flags
#    clusterCmd = 'python /storage/cylin/bin/pipeline/clusterEnhancer.py -d %s -i %s -r %s -o %s -n %s' % (chiprx_data_file,namesString,roseFolder,outputFolder,analysisName)

#    bashFile.write(clusterCmd + '\n')
#    print(clusterCmd)

#    bashFile.close()

    print('\n\n')
    print('#======================================================================')
    print('#======================XIII. MAKING GEO TABLE==========================')
    print('#======================================================================')
    print('\n\n')

    geoName = 'rasmc_chiprx'
    namesList=[]
    outputFolder = '/storage/cylin/grail/projects/rasmc_all/rasmc_geo/%s_geo/' % (geoName)

#    makeGEOCHIPRXTable(dataFile,wiggleFolder,macsFolder,namesList,geoName,outputFolder)


    print('\n\n')
    print('#======================================================================')
    print('#======================XIII. SCALING WIGGLES===========================')
    print('#======================================================================')
    print('\n\n')

    scaleTableFile = '/storage/cylin/grail/projects/rasmc_all/tables/RN6_RASMC_NEW_CHIPRX_SCALE_FACTORS.txt'
    names_list=[]
    
   # scaleWiggles(dataFile,scaleTableFile,names_list=[])




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


    crc_folder = '%scrc/rasmc_h3k27ac_0_tss_all_motifs' % (projectFolder)
    analysis_name = 'rasmc_h3k27ac_0_tss_all_motifs'
    h3k_24_list = ['RASMC_H3K27AC_PDGF_24H_NEW','RASMC_H3K27AC_PDGF_24H_REP2']
    h3k_unstim_list=['RASMC_H3K27AC_UNSTIM_NEW','RASMC_H3K27AC_UNSTIM_REP1']
#    tf_edge_brd4_delta(crc_folder,chiprx_data_file,analysis_name,h3k_24_list,h3k_unstim_list)



####################################
    # dataDict=pipeline_dfci.loadDataTable(chiprx_data_file)
    # bash_file_name='/storage/cylin/grail/projects/rasmc_all/move_chiprx_data.sh'
    # bashFile=open(bash_file_name,'w')
    # bashFile.write('#!/usr/bin/bash\n')
    # names=dataDict.keys()
    # for name in names:
    #     bam_name=dataDict[name]['bam']
    #     cmd='cp %s /storage/cylin/grail/projects/rasmc_all/bams/\n' % (bam_name)
    #     bashFile.write(cmd)


    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping chiprx bams to peaks')
    dataFile = chiprx_data_file
    dataDict=pipeline_dfci.loadDataTable(chiprx_data_file)
    names=dataDict.keys()
#    for name in names:
#        if len(dataDict[name]['enrichedMacs'])>4:
#            peak_name=dataDict[name]['enrichedMacs']
#            peak_path='%s%s' % (macsEnrichedFolder,peak_name)
#            gff_path='%s%s.gff' % (gffFolder,peak_name.split('.bed')[0])
#            utils.bedToGFF(peak_path,output=gff_path)
#            gffList=[gff_path]
#            namesL=[name]
#            pipeline_dfci.mapBamsBatch(dataFile, gffList,mappedFolder,overWrite=False,namesList=namesL,extension=0,rpm=False)


    namesL=names
    tss_gff_path='%sRN6_TSS_ALL_-300_+300.gff' % (gffFolder)
    gffList=[tss_gff_path]
    pipeline_dfci.mapBamsBatch(dataFile, gffList,mappedFolder,overWrite=False,namesList=namesL,extension=0,rpm=False)




   
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~WRITING SCALE FACTORS~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def writeScaleFactors(dataFile,namesList=[],primary_genome='hg19',rx_genome='dm6',output=''):

    '''
    creates a table of scale factors based on the rx genome read depth
    '''
    
    #first set up the output folder
    #rpm scale factor is what the rpm/bp should be MULTIPLIED by
    #mouse mapped reads give the denominator for what raw r/bp should be divided by
    outputTable = [['NAME','RN6_MAPPED_READS','DM6_MAPPED_READS','RPM_SCALE_FACTOR']]

    
    dataDict=pipeline_dfci.loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = [name for name in dataDict.keys()]
    namesList.sort()
    print('scaling the following datasets')


    for name in namesList:
        
        print('WORKING ON %s' % (name))
        bam_path = dataDict[name]['bam']
        bam = utils.Bam(bam_path)
        bam_mmr = float(bam.getTotalReads())/1000000
        scale_path = string.replace(bam_path,primary_genome,rx_genome)
        scaleBam = utils.Bam(scale_path)
        scale_mmr = float(scaleBam.getTotalReads())/1000000
        #print(bam_path)
        #print(scale_path)
        rpm_scale = bam_mmr/scale_mmr
        scale_line = [bam_mmr,scale_mmr,rpm_scale]
        scale_line = [round(x,4) for x in scale_line]
        outputTable.append([name] + scale_line)

    if len(output) == 0:
        return outputTable
    else:
        utils.unParseTable(outputTable,output,'\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~SCALING WIGGLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def scaleWiggles(dataFile,scaleTableFile,names_list=[]):

    '''
    first unzips each wiggle
    then scales each line by the rpm
    and rounds to a reasonable number (2 decimal)
    '''

    dataDict=pipeline_dfci.loadDataTable(dataFile)
    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('WCE') ==0 and name.count('INPUT') == 0]
    names_list.sort()
    print(names_list)

    print('loading scale factors')

    scale_table = utils.parseTable(scaleTableFile,'\t')
    scale_dict = {}
    for line in scale_table[1:]:
        scale_dict[line[0]] = float(line[2])
    os.chdir(wiggleFolder)
                   
#    for name in names_list:

#        print('scaling %s' % (name))
#        scale_factor = scale_dict[name]
#        wig_path_gz = '%swiggles/%s_treat_afterfiting_all.wig.gz' % (projectFolder,name)
#        wig_path = '%swiggles/%s_treat_afterfiting_all.wig' % (projectFolder,name)
        

#        wig_out = '%swiggles/%s_scaled.wig' % (projectFolder,name)
#        wig_out_final ='%swiggles/%s_scaled.wig.gz' % (projectFolder,name)
#        if utils.checkOutput(wig_out_final,0,0):
#            print('Found scaled wiggle for %s at %s' % (name,wig_out_final))
#            continue
#        cmd = 'gunzip %s' % (wig_path_gz)
#        print(cmd)
        
        
        #this should run to completion
#        os.system(cmd)

        #now open up the new wig
#        wig = open(wig_path,'r')
#        wig_scaled = open(wig_out,'w')

#        ticker = 0
#        for line in wig:

#            if ticker % 1000000 == 0:
#                print(ticker)
#            ticker+=1
#            if line[0] == 't' or line[0] == 'v':
#                wig_scaled.write(line)
#            else:
#                line = line.rstrip().split('\t')
#                line[1] = str(round(float(line[1])/scale_factor,2))
#                line_string = '\t'.join(line) + '\n'
#                wig_scaled.write(line_string)


#        wig.close()
#        wig_scaled.close()
#        cmd = 'gzip %s' % (wig_out)
#        print(cmd)
#        os.system(cmd)

#        cmd = 'gzip %s' % (wig_path)
#        print(cmd)
#        os.system(cmd)

    #now for WCE wiggles
    
    for name in names_list:

        bg_name = dataDict[name]['background']
        print('scaling %s' % (bg_name))
        scale_factor = scale_dict[bg_name]
        wig_path_gz = '%smacsFolder/%s/%s_MACS_wiggle/control/%s_control_afterfiting_all.wig.gz' % (projectFolder,name,name,name)
        wig_path = '%smacsFolder/%s/%s_MACS_wiggle/control/%s_control_afterfiting_all.wig' % (projectFolder,name,name,name)


        wig_out = '%swiggles/%s_scaled.wig' % (projectFolder,dataDict[name]['background'])
        wig_out_final ='%swiggles/%s_scaled.wig.gz' % (projectFolder,dataDict[name]['background'])
        if utils.checkOutput(wig_out_final,0,0):
            print('Found scaled wiggle for %s at %s' % (name,wig_out_final))
            continue
        cmd = 'gunzip %s' % (wig_path_gz)
        print(cmd)


        #this should run to completion
        os.system(cmd)

        #now open up the new wig
        wig = open(wig_path,'r')
        wig_scaled = open(wig_out,'w')

        ticker = 0
        for line in wig:

            if ticker % 1000000 == 0:
                print(ticker)
            ticker+=1
            if line[0] == 't' or line[0] == 'v':
                wig_scaled.write(line)
            else:
                line = line.rstrip().split('\t')
                line[1] = str(round(float(line[1])/scale_factor,2))
                line_string = '\t'.join(line) + '\n'
                wig_scaled.write(line_string)


        wig.close()
        wig_scaled.close()
        cmd = 'gzip %s' % (wig_out)
        print(cmd)
        os.system(cmd)

        cmd = 'gzip %s' % (wig_path)
        print(cmd)
        os.system(cmd)


    os.chdir(projectFolder)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING REGION GFFS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makeStitchedGFF(dataFile,set_name,names_list,write_individual=False):

    '''
    makes a stitched gff and dumps it into the gff folder
    '''

    dataDict = pipeline_dfci.loadDataTable(dataFile)
    
    loci = []
    collection_dict = {}
    gff_list = []
    for name in names_list:
        print(name)
        macsEnrichedFile = '%s%s_peaks.bed' % (macsEnrichedFolder,name)
        collection = utils.importBoundRegion(macsEnrichedFile,name)
        collection_dict[name]=collection
        loci+= collection.getLoci()
        if write_individual:
            gff_path = '%s%s_%s_-0_+0.gff' % (gffFolder,genome.upper(),name)
            gff = utils.locusCollectionToGFF(collection)
            utils.unParseTable(gff,gff_path,'\t')
            gff_list.append(gff_path)

    all_collection = utils.LocusCollection(loci,50)
    stitched_collection = all_collection.stitchCollection()

    gff = utils.locusCollectionToGFF(stitched_collection)

    out_path = '%s%s_%s_UNION_-0_+0.gff' % (gffFolder,genome.upper(),set_name)
    print(out_path)
    utils.unParseTable(gff,out_path,'\t')
    gff_list.append(out_path)

    #now get the intersect gff
    print('getting intersection gff')
    stitched_loci = stitched_collection.getLoci()
    intersect_loci = []
    ticker = 0
    for locus in stitched_loci:
        if ticker%1000==0:
            print(ticker)
        ticker+=1
        overlap = True
        for name in names_list:
            if len(collection_dict[name].getOverlap(locus,'both')) == 0:
                overlap = False

        if overlap == True:
            intersect_loci.append(locus)

            
    print('identified %s stitched loci' % (len(stitched_loci)))
    print('identified %s intersect loci' % (len(intersect_loci)))

    intersect_collection = utils.LocusCollection(intersect_loci,50)

    intersect_gff = utils.locusCollectionToGFF(intersect_collection)

    intersect_path = '%s%s_%s_INTERSECT_-0_+0.gff' % (gffFolder,genome.upper(),set_name)
    print(intersect_path)
    utils.unParseTable(intersect_gff,intersect_path,'\t')
    gff_list.append(intersect_path)
    return gff_list



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~MAKING ACTIVE GENE LIST~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def makeActiveTable(h3k_mapped='',expr_table=''):
    #Takes in path of h3k27ac mapped enriched file for all h3k27ac bams at the TSS +- 1000 bp
    #and the path to the rna expression fpkm means table
    #Any gene found in both is included in the list of active genes

    h3k27ac_mapped_Table = utils.parseTable(h3k_mapped,'\t')
    fpkm_expr_table = utils.parseTable(expr_table,'\t')
    activeList = []

    for x in h3k27ac_mapped_Table:
        for y in fpkm_expr_table:
            if y[0] == x[2]:
                if float(y[1]) >= 10 or float(y[2]) >= 10 or float(y[4]) >=10:
                    new_line = [x[1],x[2]]
                    activeList.append(new_line)

    out_path = '%stables/activeListTable_expr_filter.txt' % (projectFolder)
    activeListTable = utils.unParseTable(activeList,out_path,'\t')
    print('active table located at %s') % (out_path)


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




#============================================================================
#================MAKE CHIPRX GEO TABLE=======================================
#============================================================================

def makeGEOCHIPRXTable(dataFile,wiggleFolder,macsFolder,namesList,geoName,outputFolder =''):

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
    bashFile.write('\n\n')


    #write the untar command
    for name in namesList:
        fastqFile = dataDict[name]['fastq']
        if len(fastqFile) == 0:
            print("NO FASTQ FILE FOR %s") % (name)
            continue
        if fastqFile.count('tar.gz') > 0: #for files generated by whitehead that have tar header #####RACHEL READ HERE
            tarCmd = 'tar --strip-components 5 --to-stdout -xzvf %s | gzip > %s.fastq.gz\n' % (fastqFile,name)
        else:
            tarCmd = 'cp %s %s.fastq.gz\n' % (fastqFile,name)
        bashFile.write(tarCmd)

    bashFile.write('\n\n\n')


    #write the wiggle cp command
    for name in namesList:
        if name.count('WCE') == 1 or name.count('INPUT') == 1:
            refName = backgroundDict[name]
            controlWiggleFile = '%s%s/%s_MACS_wiggle/control/%s_control_afterfiting_all.wig.gz' % (macsFolder,refName,refName,refName)
            wigCmd = "cp '%s' %s.wig.gz\n" % (controlWiggleFile,name)
            #wigCmd = "cp '%swceWiggles/%s_control_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,refName,name)
        else:
            wigCmd = "cp '%s%s_treat_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,name,name)
        bashFile.write(wigCmd)


    for name in namesList:
        md5Cmd = 'md5sum %s.chiprx.scaled.bedgraph.gz >> md5sum.txt\n' % (name)
        bashFile.write(md5Cmd)

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
#~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF CLUSTERS~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def tf_edge_brd4_delta_out(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list,output=''):

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
    print('making log2 24 vs unstim signal table at edges')
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

    print('unstim median H3k27ac signal')
    print(o_median)
    print('24H median H3k27ac signal')
    print(y_median)

    #now that we have the median, we can take edges where at least 1 edge is above the median
    #and both are above zero and generate a new table w/ the fold change

    signal_filtered_path = string.replace(signal_path,'.txt','_filtered.txt')
    if utils.checkOutput(signal_filtered_path,0,0):
        print('Found filtered signal table for edges at %s' % (signal_filtered_path))
        signal_table_filtered = utils.parseTable(signal_filtered_path,'\t')
    else:

        signal_table_filtered = [signal_table[0]+['UNSTIM_MEAN','24H_MEAN','24_vs_Unstim_LOG2']]
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
    print('making log2 y vs o signal table at edges')
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




#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
