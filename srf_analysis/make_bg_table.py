import sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils
import os

motifBedDir=sys.argv[1]

list_files=os.listdir(motifBedDir)
allLoci=[]
bg_files=[]
for f in list_files:
    fEnd=f.split('.')[-1]
    if fEnd == 'meme':
        motif_beds.append(f)

all_bg=[]

for bg in bg_files:
	bg_tab=utils.parseTable(bg,' ')
	print(bg_tab[0])
	print(bg_tab[8])
	new_line=[]
	if len(all_bg) == 0:
		new_line=['TF']
		for line in bg_tab[8:]:
			new_line.extend(line[0])
		all_bg.append(new_line)
	else:
		tf_name=(bg.split('/')[-1]).split('_')[0]
		new_line=[tf_name]
		for line in bg_tab[8:]:
			new_line.extend(line[1])
		all_bg.append(new_line)

utils.unParseTable(all_bg,'/storage/cylin/grail/projects/rasmc_all/tables/all_motif_bgs.txt','\t')