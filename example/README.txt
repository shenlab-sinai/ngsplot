This folder contains two toy examples for you to get a sense of ngs.plot workflow. 
User can also use these examples to check the completion of installation.

========Example for ChIP-seq data:========
chipseq_mm9_ch19.bam is the alignment file (mm9, chr19 only) for histone mark 
H3K4me3 in BAM format, derived from mouse neurons. Here we demonstrate how to 
draw a TSS plot in two simple steps:

1. Calculate the coverage file based on the alignment file.
cov.calc.r -G mm9 -R chipseq_mm9_chr19.bam -X bam -C chipseq_mm9_chr19.RData

2. Plot a figure at TSS.
ngs.plot.r -R tss -C chipseq_mm9_chr19.RData -O chipseq_mm9_chr19 -T ChIP-seq:MM9:chr19 -SE 0

It takes ~20 sec to calculate coverage and ~2 min to draw the plot. Users are 
encouraged to experiment with various visualization options to feel the 
convenience and power of ngs.plot.


========Example for RNA-seq data:========
rnaseq_rn4_chr20.bam is an RNA-seq alignment file (rn4, chr20 only). We are 
going to draw a genebody plot in two simple steps:

1. Calculate the coverage file, same as above.
cov.calc.r -G rn4 -R rnaseq_rn4_chr20.bam -X bam -C rnaseq_rn4_chr20.RData

2. Plot a figure on genebody.
ngs.plot.r -R genebody -C rnaseq_rn4_chr20.RData -F rnaseq -O rnaseq_rn4_chr20 -T RNA-seq:RN4:chr20 -SE 0

It takes ~15 sec to calculate coverage and ~1 min to draw the plot. Notice that 
we used "-F rnaseq" here to concatenate coverage from exons. This is necessary 
for RNA-seq because introns barely contain any coverage information. If you do 
not use this option, the figure will simply be ugly.


========Example configuration and gene list file:========

ngs.plot can also be used to generate one figure with multiple samples using a 
configuration file. An example is given in "config.cov.example.txt".

Gene list can be used to subset genomic coverage. A list can be either used in
command line or supplied in a configuration file. The format is simply one
gene/transcript symbol/ID per line. Both RefSeq and Ensembl annotations are
accepted. You can even mix heterogeneous symbols/IDs. An example is given in 
"example.gene.list.txt".



For more details of the usage, please visit project website:
http://code.google.com/p/ngsplot/








