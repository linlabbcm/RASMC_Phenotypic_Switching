import sys
sys.path.append('/storage/cylin/bin/pipeline')
import utils

print "Running: ",str(sys.argv)

bed_file_path=sys.argv[1]

bed=utils.parseTable(bed_file_path,'\t')

new_bed=[]
i=0
for line in bed:
	print i
	i=i+1
	chrom=line[0]
	start=int(line[1])
	stop=int(line[2])
	reg_id=line[3]
	mod=(stop-start)%50
	if mod%2 == 0:
		n_start=start-(mod/2)
		n_stop=stop+(mod/2)
	else:
		mod=mod-1
		start=start-((mod/2)+1)
		stop=stop+(mod/2)
	count=n_start
	while count > stop:
		new_line=[chrom,count,count+50,reg_id]
		new_bed.append(new_line)
		count=count+50

a = bed_file_path.split(".")[0]
new_file_path='{}{}'.format(a,'_50bp_bins.bed')
utils.unParseTable(new_bed,new_file_path,'\t')