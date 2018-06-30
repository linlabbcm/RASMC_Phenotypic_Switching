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
        bg_files.append(f)

all_bg=[]

for bg in bg_files:
        filename='%s%s' % (motifBedDir,bg)
	bg_tab=utils.parseTable(filename,' ')
	print(bg_tab[0])
	print(bg_tab[8])
	new_line=[]
	if len(all_bg) == 0:
		new_line=['TF']
		for line in bg_tab[8:]:
			new_line.append(line[0])
		all_bg.append(new_line)

        tf_name=(bg.split('/')[-1]).split('_')[0]
        new_line=[tf_name]
        for line in bg_tab[8:]:
                new_line.append(line[1])
        all_bg.append(new_line)

utils.unParseTable(all_bg,'/storage/cylin/grail/projects/rasmc_all/tables/all_motif_bgs_all_TFs.txt','\t')
