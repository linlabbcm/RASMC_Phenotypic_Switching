import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils
import re


gff_path = '/storage/cylin/grail/projects/rasmc_all/gff/rasmc_h3k27ac_0_tss_all_subpeak.gff'
genome_directory='/storage/cylin/grail/genomes/Rattus_norvegicus/UCSC/rn6/Sequence/Chromosomes/'

genome = 'RN6'

print('gffToFasta Tool running on ' + gff_path + ' for ' + genome)
fasta = utils.gffToFasta(genome,genome_directory,gff_path,UCSC=True,useID=False)


print('Creating density table')
table=[]
header=['DENSITY','POSITIONS','POS_COUNT','SUBPEAK_LENGTH']
table.append(header)

#CArG box motif
seq='CC[AT]{6}GG'

table_path='/storage/cylin/grail/projects/rasmc_all/motif_density/CArG_box_seq_density_from_fasta_full_length_no_slash.txt'




for i in range(0,len(fasta),2):
    positions=[]
    line=fasta[i+1]
    forward_count=re.findall(seq,line)
    f_pos = re.finditer(seq,line)
    for j in f_pos:
        positions.append(str(j.span()[0]))
        print(positions)
    
    seq_density = len(positions)*10/len(line)
    pos_string = ','.join(positions)
    new_line=[seq_density,pos_string,len(line)]
    table.append(new_line)
#    print(seq_density)


print('Printing out table '+table_path)

utils.unParseTable(table,table_path,'\t')

scr_table=[]
scr_table.append(header)
for i in range(0,len(fasta),2):
    line=fasta[i+1]
    scramble_count=0
    positions=[]
    for j in range(len(line)-10):
        k_mer = line[j:j+10]
        if k_mer.count('C')==2 and k_mer.count('G')==2:
            if (k_mer[0:2]=='CC' and k_mer[-2:]=='GG'):
                continue
            else:
                scramble_count+=1
                positions.append(str(j))
    pos_string = ','.join(positions)
    new_line=[float(((len(positions))*10)/len(line)),pos_string,len(positions),len(line)]
    scr_table.append(new_line)

scr_table_path='/storage/cylin/grail/projects/rasmc_all/motif_density/CArG_box_scramble_density_from_fasta_full_length_scr.txt'

print('Printing out table '+scr_table_path)

utils.unParseTable(scr_table,scr_table_path,'\t')




