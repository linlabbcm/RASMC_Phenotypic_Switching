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

#rna
rna_data_file = '%sdata_tables/RASMC_RNA_DATA_TABLE.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('rna analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I, LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(rna_data_file)

    
    print('\n\n')
    print('#==========================================================================')
    print('#=======================II, ALIGNING WITH HISAT2===========================')
    print('#==========================================================================')
    print('\n\n')

    
    #pipeline_dfci.mapHisat(dataFile,namesList=[],useSRA=False,pCount=16,Launch=True)



    print('\n\n')
    print('#==========================================================================')
    print('#=======================III, RUNNING RNA-SEQ ANALYSIS======================')
    print('#==========================================================================')
    print('\n\n')

    #analysisName = 'rasmc_rna'

    #gtfFile = '/storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf'
    #cufflinksFolder = utils.formatFolder('%scufflinks' % (projectFolder),True)
    #bashFileName = '%s%s_cufflinks.sh' % (cufflinksFolder,analysisName)

    #groupList = [['RASMC_RNA_0H_A','RASMC_RNA_0H_B'],['RASMC_RNA_PDGF_2H_B','RASMC_RNA_PDGF_2H_C','RASMC_RNA_PDGF_2H_D'],['RASMC_RNA_PDGF_JQ1_2H_E','RASMC_RNA_PDGF_JQ1_2H_G','RASMC_RNA_PDGF_JQ1_2H_H'],['RASMC_RNA_PDGF_24H_A','RASMC_RNA_PDGF_24H_B','RASMC_RNA_PDGF_24H_D'],['RASMC_RNA_PDGF_JQ1_24H_E','RASMC_RNA_PDGF_JQ1_24H_F','RASMC_RNA_PDGF_JQ1_24H_H']]

    #print(groupList)
    #pipeline_dfci.makeCuffTableSlurm(rna_data_file,analysisName,gtfFile,cufflinksFolder,groupList,bashFileName)

    # #flag useERCC to true


    print('\n\n')
    print('#==========================================================================')
    print('#=======================IV, MAKE GEO TABLE=================================')
    print('#==========================================================================')
    print('\n\n')



    namesList = []
    geoName = 'rasmc_rna'

    outputFolder = '/storage/cylin/grail/projects/rasmc_all/rasmc_geo/%s_geo/' % (geoName)
    makeGEORNATable(rna_data_file,namesList,geoName,outputFolder)





#======================================================================================================
#==================================FUNCTIONS FOR ANALYSIS==============================================
#======================================================================================================

#==========================================================================
#===========================MAKE GEO RNA TABLE=============================
#==========================================================================


def makeGEORNATable(dataFile,namesList,geoName,outputFolder =''):

    '''
    makes a geo table and a bash script to format everything
    '''
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    if len(namesList) == 0:
        namesList = dataDict.keys()

    #set up bash script and output folder
    outputFolder = pipeline_dfci.formatFolder(outputFolder,True)
    bashFileName = '%s%s_bash.sh' % (outputFolder,geoName)

    bashFile = open(bashFileName,'w')


    geoTable = [['SAMPLE_NAME','TITLE','CELL_TYPE','RAW_FILE','BARCODE']]


    namesList.sort()

    for name in namesList:

        sampleName = dataDict[name]['uniqueID']
        title = name
        cell_type = name.split('_')[0]
        raw_file = "%s.fastq.gz" % (name)

        fastqFile = dataDict[name]['fastq']
        uniqueID = dataDict[name]['uniqueID']
        try:
            barcode = pipeline_dfci.getTONYInfo(uniqueID,38)
        except IndexError:
            barcode = ''

        newLine = [sampleName,title,cell_type,raw_file,barcode]
        geoTable.append(newLine)


    utils.unParseTable(geoTable,"%s%s_meta.xls" % (outputFolder,geoName),'\t')

    #now make the folder to hold everything and the relevant bash script
    if len(outputFolder) == 0:
        outputFolder ='./%s/' % (geoName)

    else:
        outputFolder = outputFolder + geoName + '/'

    pipeline_dfci.formatFolder(outputFolder,True)



    #now make the bash file
    bashFile.write('#!/usr/bin/bash\n')
    bashFile.write('cd %s\n' %(outputFolder))
    bashFile.write('\n\n')


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



    #write the md5sums for the wiggles
    bashFile.write('\n\n\n')
    bashFile.write("echo '' > md5sum.txt\n")

    #write md5sums for the fastqs
    for name in namesList:
        md5Cmd = 'md5sum %s.fastq.gz >> md5sum.txt\n' % (name)
        bashFile.write(md5Cmd)

    #the big tar command
    tarCmd = '#tar -cvzf %s.tar.gz %s\n' % (geoName,outputFolder)
    bashFile.write(tarCmd)
    bashFile.close()




#==========================================================================
#==================================THE END=================================
#==========================================================================


if __name__=="__main__":
    main()







