import os,sys
sys.path.append('/storage/cylin/bin/pipeline/')
import utils


motif_bed_dir = '/storage/cylin/grail/projects/rasmc_all/beds/srf_motif_analysis/'
motif_beds = os.listdir(motif_bed_dir)

allLoci = []

for bed in motif_beds:
	TF_name = bed.split('_')[0]
	collection = utils.importBoundRegion('%s%s' %(motif_bed_dir,bed),TF_name)

	allLoci += collection.getLoci()



giant_collection = utils.LocusCollection(allLoci,50)

stitched_collection = giant_collection.stitchCollection(stitchWindow=50)

new_bed = utils.locusCollectionToBed(stitched_collection)

utils.unParseTable(new_bed,'%s50_bp_stitched_srf_motif_analysis_bed.bed' % (motif_bed_dir),'\t')



