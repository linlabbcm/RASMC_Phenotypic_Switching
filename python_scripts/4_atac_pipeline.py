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
macsFolder = '%s/riesling/peaks/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%s/riesling/wiggles/' % (projectFolder)
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

#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#atac-seq
#atac_data_file = '%sdata_tables/RASMC_ATAC_TABLE.txt' % (projectFolder)
atac_data_file = '/storage/cylin/grail/projects/rasmc/RASMC_ATAC_TABLE.txt'




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
    pipeline_dfci.summary(atac_data_file)

    #assumes macs has already been run and formatted
#    run_macs(chip_data_file)

#    sys.exit()



    print('\n\n')
    print('#======================================================================')
    print('#======================MAKING GEO TABLES===============================')
    print('#======================================================================')
    print('\n\n')


    geoName='rasmc_atac'
    outputFolder = '/storage/cylin/grail/projects/rasmc_all/rasmc_geo/%s_geo/' % (geoName)
    namesList=[]

#    makeGEOTable(atac_data_file,wiggleFolder,macsFolder,namesList,geoName,outputFolder)


    #==========================================================================
    #====================MAP BAMS BATCH========================================
    #==========================================================================
    print('Mapping chiprx bams to peaks')
    dataFile = atac_data_file
    dataDict=pipeline_dfci.loadDataTable(atac_data_file)
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
    pipeline_dfci.mapBamsBatch(dataFile, gffList,mappedFolder,overWrite=False,namesList=namesL,extension=0,rpm=True)





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
        tarCmd = 'cp %s %s.fastq.gz\n' % (fastqFile,name)
        bashFile.write(tarCmd)

    bashFile.write('\n\n\n')
    #write the wiggle cp command
    for name in namesList:
        if name.count('WCE') == 1 or name.count('INPUT') == 1 :
            refName = backgroundDict[name]
            controlWiggleFile = '%s%s/%s_MACS_wiggle/control/%s_control_afterfiting_all.wig.gz' % (macsFolder,refName,refName,refName)
            wigCmd = "cp '%s' %s.wig.gz\n" % (controlWiggleFile,name)
            #wigCmd = "cp '%swceWiggles/%s_control_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,refName,name)
        else:
            wigCmd = "cp '%s%s_treat_afterfiting_all.wig.gz' %s.wig.gz\n" % (wiggleFolder,name,name)
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


#==========================================================================
#==================================THE END=================================
#==========================================================================


if __name__=="__main__":
    main()

