This folder contains the example files of ngs.plot. Users could use these examples to check the installation, and have an overview of ngs.plot workflow.

Example for ChIP-seq data:
chipseq_mm9_ch19.bam is an H3K4me3 alignment file of chr19, mm9. Here we demostrate how to plot a TSS plot of it.

# To calculate the coverage file based on given alignment file.
cov.calc.r -G mm9 -R chipseq_mm9_chr19.bam -X bam -C chipseq_mm9_chr19.RData
# To plot the figure at TSS.
ngs.plot.r -R tss -C chipseq_mm9_chr19.RData -O chipseq_mm9_chr19 -T ChIP-seq:MM9:chr19 -SE 0

Example for RNA-seq data:
rnaseq_rn4_chr20.bam is an RNA-seq data of chr20, rn4. We could draw a plot of genebody based on the alignment.

# To calculate the coverage file, same step as above.
cov.calc.r -G rn4 -R rnaseq_rn4_chr20.bam -X bam -C rnaseq_rn4_chr20.RData
# To plot the figure of genebody.
ngs.plot.r -R genebody -C rnaseq_rn4_chr20.RData -F rnaseq -O rnaseq_rn4_chr20 -T RNA-seq:RN4:chr20 -SE 0

ngs.plot could also be used to generate one plot with multiple datasets by configuration file, or based on given gene lists. 
Example of gene list is example.gene.list.txt.
The configure file for multiple datasets is config.cov.example.txt.

For more details of the usage, visit project website please:
http://code.google.com/p/ngsplot/
