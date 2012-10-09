if [ ! -f ./merged.bed ]
then
gtf2tree.pl   ~/reference/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf |sort|uniq >  merged.bed
fi


get_difflist.pl  isoform_exp.diff  Control:Treat  up.isoform_exp.diff   down.isoform_exp.diff

parse_diff.pl  splicing.diff  splicing.significant.diff

for i in variant  altBoth  altDonor  altAcceptor; do alter2bed.pl merged.bed  splicing.significant.diff  up.isoform_exp.diff  ${i}  up.alternative_splicing.${i}.significant.bed.tmp; awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}'  up.alternative_splicing.${i}.significant.bed.tmp |sort|uniq > up.alternative_splicing.${i}.significant.bed; rm up.alternative_splicing.${i}.significant.bed.tmp ;done

for i in variant  altBoth  altDonor  altAcceptor; do alter2bed.pl merged.bed  splicing.significant.diff down.isoform_exp.diff  ${i}  down.alternative_splicing.${i}.significant.bed.tmp; awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}'  down.alternative_splicing.${i}.significant.bed.tmp |sort|uniq > down.alternative_splicing.${i}.significant.bed; rm  down.alternative_splicing.${i}.significant.bed.tmp;done


### rm the tmp files
rm  merged.bed   up.isoform_exp.diff   down.isoform_exp.diff   splicing.significant.diff  
