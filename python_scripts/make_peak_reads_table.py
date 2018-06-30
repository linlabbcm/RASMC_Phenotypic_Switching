import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils
import re
import os
import subprocess
import numpy

mapped_peak_dir='/storage/cylin/grail/projects/rasmc_all/mappedFolder/mapped_peaks/'

list_files=os.listdir(mapped_peak_dir)

reads_table=[]

for f in list_files:
	f_path='%s%s' % (mapped_peak_dir,f)
	f_table=utils.parseTable(f_path,'\t')
	peak_rpm=0
	for line in f_table[1:]:
		start=float((line[1].split(':')[1]).split('-')[0])
		stop=float((line[1].split(':')[1]).split('-')[1])
		length=stop-start
		AUC=float(line[2])
		peak_rpm=peak_rpm+(AUC*(length/2))
	new_line=[f,peak_rpm]
	reads_table.append(new_line)

utils.unParseTable(reads_table, '/storage/cylin/grail/projects/rasmc_all/ATAC_peak_reads_table.txt','\t')


mapped_peak_dir='/storage/cylin/grail/projects/rasmc_all/mappedFolder/tss_peaks/'

list_files=os.listdir(mapped_peak_dir)

reads_table=[]

for f in list_files:
        f_path='%s%s' % (mapped_peak_dir,f)
        f_table=utils.parseTable(f_path,'\t')
        peak_rpm=0
        for line in f_table[1:]:
                start=float((line[1].split(':')[1]).split('-')[0])
                stop=float((line[1].split(':')[1]).split('-')[1])
                length=stop-start
                AUC=float(line[2])
                peak_rpm=peak_rpm+(AUC*(length/2))
        new_line=[f,peak_rpm]
        reads_table.append(new_line)

utils.unParseTable(reads_table, '/storage/cylin/grail/projects/rasmc_all/ATAC_TSS_reads_table.txt','\t')

mapped_peak_dir='/storage/cylin/grail/projects/rasmc_all/mappedFolder/tss_peaks/'

list_files=os.listdir(mapped_peak_dir)

reads_table=[]

for f in list_files:
        f_path='%s%s' % (mapped_peak_dir,f)
        f_table=utils.parseTable(f_path,'\t')
        peak_rpm=0
        for line in f_table[1:]:
                start=float((line[1].split(':')[1]).split('-')[0])
                stop=float((line[1].split(':')[1]).split('-')[1])
                length=stop-start
                AUC=float(line[2])
                peak_rpm=peak_rpm+(AUC*length)
        new_line=[f,peak_rpm]
        reads_table.append(new_line)

utils.unParseTable(reads_table, '/storage/cylin/grail/projects/rasmc_all/TSS_reads_table.txt','\t')



