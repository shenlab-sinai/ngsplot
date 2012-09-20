if [ ! -f ./merged.bed ]
then
gtf2tree.pl   ~/reference/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf   merged.bed
fi


get_difflist.pl  isoform_exp.diff  Control:Treat  up.isoform_exp.diff   down.isoform_exp.diff

parse_diff.pl  splicing.diff  splicing.significant.diff

for i in variant  altBoth  altDonor  altAcceptor; do alter2bed.pl merged.bed  splicing.significant.diff  up.isoform_exp.diff  ${i}  up.alternative_splicing.${i}.significant.bed;done

for i in variant  altBoth  altDonor  altAcceptor; do alter2bed.pl merged.bed  splicing.significant.diff down.isoform_exp.diff  ${i}  down.alternative_splicing.${i}.significant.bed;done



