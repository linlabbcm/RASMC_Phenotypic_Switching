import os,sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils



motif_bed_dir = '/storage/cylin/grail/projects/rasmc_all/crc/rasmc_h3k27ac_0_tss/motif_beds/'
motif_beds = os.listdir(motif_bed_dir)

bashFileName = %ssort.sh % (motif_bed_dir)
bashFile = open(bashFileName,'w')

bashFile.write('#/usr/bin/bash\n')

for bed in motif_beds:
	TF_name = bed.split('_')[0]

	bashFile.write('sort -k 1,1 -k2,2n %s%s_motifs.bed > %s%s_motifs_sorted.bed \n' % (motif_bed_dir,TF_name,motif_bed_dir,TF_name))

bashFile.close()	