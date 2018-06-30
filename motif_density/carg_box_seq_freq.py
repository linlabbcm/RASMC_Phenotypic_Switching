import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils
import re



#def motifFrequency(analysis_name,motifBedDir,genome,seq,output,tf_list_path=''):
    
analysis_name=sys.argv[1]
motifBedDir=sys.argv[2]
genome=sys.argv[3]
seq=sys.argv[4]
output=sys.argv[5]
if len(sys.argv) > 6:
    tf_list_path=sys.argv[6]
else:
    tf_list_path=''

genome_dir_dict={'RN6':'/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Chromosomes/'}

list_files=os.listdir(motifBedDir)
allLoci=[]
motif_beds=[]
for f in list_files:
    fEnd=f.split('.')[-1]
    if fEnd == 'bed':
        motif_beds.append(f)

print(len(motif_beds))

temp_dir = '%stmp/' % (motifBedDir)
figures_dir = '%sfigures/' % (motifBedDir)
tables_dir = '%stables/' % (motifBedDir)

#making folders
folderList = [temp_dir,figures_dir,tables_dir]

for folder in folderList:
    formatFolder(folder,True)

print(tf_list_path)
print(len(tf_list_path))

tf_list=[]
if len(tf_list_path)>0:
    tf_table=utils.parseTable(tf_list_path,'\t')
    for tf in tf_table:
        print(tf[0])
        tf_list.append(tf[0].upper())


    print(tf_list)

    tf_beds=[]
    for bed in motif_beds:
        tf_name=bed.split('_')[0]
        if tf_name in tf_list:
            tf_beds.append(bed)
            print(tf_beds)

    motif_beds=tf_beds


freq_table=[]
header=['TF_NAME','FREQ_SUM','NUM_REGIONS','MOTIF/KB']
#remove 'track' line from crc bed files and write to tmp file for Rscript use
for bed in motif_beds:
    bed_path = '%s%s' % (motifBedDir,bed)
    print(bed)
    bed_table = utils.parseTable(bed_path,'\t')
    new_table=[]
    for line in bed_table[1:]:
        chrom=line[0]
        start=int(line[1])-50
        stop=int(line[2])+50
        edge=line[3]
        strand=line[4]
        new_line=[chrom,start,stop,edge,strand]
        new_table.append(new_line)
    tmp_path = '%s%s' % (temp_dir,bed)
    if len(new_table) > 0:
        utils.unParseTable(new_table,tmp_path,'\t')
        gff=bed.split('.')[0]+'.gff'
        print(gff)
        gff_path='%s%s' % (temp_dir,gff)
        bedToGFF(tmp_path,gff_path)

    



        genome_directory=genome_dir_dict[genome]



        print('gffToFasta Tool running on ' + gff_path + ' for ' + genome)
        fasta = utils.gffToFasta(genome,genome_directory,gff_path,UCSC=True,useID=False)


        print('Creating density table')
        table=[]
        header=['DENSITY','POSITIONS','SUBPEAK_LENGTH']
        table.append(header)

        #CArG box motif 'CC[A\T]{6}GG'

        table_path=output




        for i in range(0,len(fasta),2):
            positions=[]
            line=fasta[i+1]
            forward_count=re.findall(seq,line)
            f_pos = re.finditer(seq,line)
            for j in f_pos:
                positions.append(str(j.span()[0]))
#               print(positions)
    
            seq_density = float(((len(positions))*10)/len(line))
            pos_string = ','.join(positions)
            new_line=[seq_density,pos_string,len(line)]
            table.append(new_line)
#           print(seq_density)
        tab_len=len(table[1:])
        tf_name=bed.split('_')[0]
        print(tf_name)
        freq_sum=0
        for line in table[1:]:
            freq_sum=freq_sum+float(line[0])
        freq_tab_line=[tf_name,freq_sum,tab_len,(freq_sum/tab_len)*1000]
        freq_table.append(freq_tab_line)




print('Printing out table '+table_path)

utils.unParseTable(freq_table,table_path,'\t')




