#!/usr/bin/bash
cd /storage/cylin/grail/projects/rasmc_all/cufflinks/

echo 'making cuffquant folders'
mkdir RASMC_RNA_0H_A
mkdir RASMC_RNA_0H_B
mkdir RASMC_RNA_PDGF_2H_B
mkdir RASMC_RNA_PDGF_2H_C
mkdir RASMC_RNA_PDGF_2H_D
mkdir RASMC_RNA_PDGF_JQ1_2H_E
mkdir RASMC_RNA_PDGF_JQ1_2H_G
mkdir RASMC_RNA_PDGF_JQ1_2H_H
mkdir RASMC_RNA_PDGF_24H_A
mkdir RASMC_RNA_PDGF_24H_B
mkdir RASMC_RNA_PDGF_24H_D
mkdir RASMC_RNA_PDGF_JQ1_24H_E
mkdir RASMC_RNA_PDGF_JQ1_24H_F
mkdir RASMC_RNA_PDGF_JQ1_24H_H

# echo 'calling cuffquant'
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_0H_A/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1212.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_0H_B/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1202.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_2H_B/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1213.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_2H_C/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1203.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_2H_D/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1208.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_2H_E/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1214.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_2H_G/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1204.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_2H_H/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1209.RN6_ERCC.sorted.bam --library-type fr-firststrand 
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_24H_A/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1215.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_24H_B/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1205.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_24H_D/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1210.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_24H_E/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1216.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_24H_F/ /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1206.RN6_ERCC.sorted.bam --library-type fr-firststrand &
# cuffquant -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_24H_H/ /storage/cylin/grail/genomets/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/bam/rn6_ercc/20140328_1211.RN6_ERCC.sorted.bam --library-type fr-firststrand &

# echo 'running cuffnorm command'
# cuffnorm -p 4 -o /storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_all_rna_cuffnorm/ -L RASMC_RNA_0H,RASMC_RNA_PDGF_2H,RASMC_RNA_PDGF_JQ1_2H,RASMC_RNA_PDGF_24H,RASMC_RNA_PDGF_JQ1_24H /storage/cylin/grail/genomes/ERCC_Technical_Data/rn6_ercc.gtf /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_0H_A/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_0H_B/abundances.cxb /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_2H_B/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_2H_C/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_2H_D/abundances.cxb /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_2H_E/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_2H_G/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_2H_H/abundances.cxb /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_24H_A/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_24H_B/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_24H_D/abundances.cxb /storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_24H_E/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_24H_F/abundances.cxb,/storage/cylin/grail/projects/rasmc_all/cufflinks/RASMC_RNA_PDGF_JQ1_24H_H/abundances.cxb

R --no-save /storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_all_rna_cuffnorm/genes.fpkm_table /storage/cylin/grail/projects/rasmc_all/cufflinks/rasmc_all_rna_cuffnorm/output/ rasmc_all_rna RASMC_RNA_0H,RASMC_RNA_PDGF_2H,RASMC_RNA_PDGF_JQ1_2H,RASMC_RNA_PDGF_24H,RASMC_RNA_PDGF_JQ1_24H TRUE < /storage/cylin/bin/pipeline/normalizeRNASeq.R
