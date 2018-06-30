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
header=['DENSITY','POSITIONS','SUBPEAK_LENGTH']
table.append(header)

#Nr1d1 node motifs
#seq='AGGTCA'
#rev_seq='TCCAGT'
#table_path='/storage/cylin/grail/projects/rasmc_all/motif_density/seq_density_from_fasta_full_length.txt'


#Jun/Fos motifs
seq='TGACTCA'
rev_seq='ACTGAGT'
table_path='/storage/cylin/grail/projects/rasmc_all/tables/seq_density_from_fasta_jun_fos.txt'

for i in range(0,len(fasta),2):
    positions=[]
    line=fasta[i+1]
    forward_count=re.findall(seq,line)
    reverse_count=re.findall(rev_seq,line)
    f_pos = re.finditer(seq,line)
    r_pos = re.finditer(rev_seq,line)
    for j in f_pos:
        positions.append(str(j.span()[0]))
#        print(positions)
    for k in r_pos:
        positions.append(str(k.span()[1]))       
#        print(positions)
    seq_density = float((len(forward_count))+(len(reverse_count)))/len(line)
    pos_string = ','.join(positions)
    new_line=[seq_density,pos_string,len(line)]
    table.append(new_line)
#    print(seq_density)


print('Printing out table '+table_path)

utils.unParseTable(table,table_path,'\t')




