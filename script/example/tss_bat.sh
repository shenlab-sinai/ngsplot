if [ ! -f ./merged.bed ]
then
gtf2tree.pl  ~/reference/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf  merged.bed
fi

get_difflist.pl  tss_group_exp.diff   Control:Treat  up.tss_group_exp.diff   down.tss_group_exp.diff

parse_diff.pl  promoters.diff  promoters.significant.diff

combine_diff.pl promoters.significant.diff  up.tss_group_exp.diff  up.tss_group_exp.both.significant.diff

combine_diff.pl promoters.significant.diff  down.tss_group_exp.diff  down.tss_group_exp.both.significant.diff

difflist2bed.pl  merged.bed  up.tss_group_exp.both.significant.diff   up.tss_group_exp.both.significant.bed 

difflist2bed.pl  merged.bed  down.tss_group_exp.both.significant.diff  down.tss_group_exp.both.significant.bed

coordinat.pl up.tss_group_exp.both.significant.bed  exon:promoter  up.tss.significant.bed

coordinat.pl down.tss_group_exp.both.significant.bed  exon:promoter  down.tss.significant.bed


