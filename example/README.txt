This folder contains four toy samples for you to get a sense of the ngs.plot 
workflow. User can also use these examples to check the completion of 
installation.

Due to Google code's 200M limitation, we have to separate the main program and 
the bam example files. If you want to run through the following examples, un-
archive and merge the "example" and "example.bam" folders.


========Example for ChIP-seq data:========

hesc.H3k27me3.1M.bam and hesc.H3k4me3.1M.bam are alignment files for two 
histone marks, H3k27me3 and H3k4me3, from the ENCODE project. hesc.Input.500K.
bam is an alignment file for DNA input. For demonstration purposes, I only 
kept 1M or 500K alignments for each file.

1. First, I show a simple example of doing plot for one histone mark around 
   TSS, this can be done easily using command:

   ngs.plot.r -G hg19 -R tss -C hesc.H3k4me3.1M.bam -O k4.test

   ngs.plot.r will generate three files: k4.test.avgprof.pdf, k4.test.heatmap.
   pdf and k4.test.zip. The first two files are plots by default settings. The 
   zip file contains data which can be used for re-plotting.

2. Now, I give an example of drawing plots for two marks together. You'll need 
   to create a configuration file, such as:

   hesc.H3k4me3.1M.bam<TAB>-1<TAB>"H3K4me3"
   hesc.H3k27me3.1M.bam<TAB>-1<TAB>"H3K27me3"

   You need to press the tab key on your keyboard to replace <TAB> for your 
   real-life configuration file. This tells ngs.plot to draw a plot for the 
   two marks on the whole genome. "-1" means all genes in the database. 
   However, you can create your own gene lists (in text files) and input here 
   to subset the drawing regions. Assume your configuration file is named 
   "config.example.txt", the command is like:

   ngs.plot.r -G hg19 -R tss -C config.example.txt -O encode1M.k4k27

3. Each time you perform a study with ngs.plot, it will generate a zip file.
   This file is very helpful because you can use it to re-create images with 
   different parameter values. It is also useful when you do heavy-duty 
   computations on a remote server and want to download the results to your 
   local computer for visualization. Again, I use the H3k4me3 data as an 
   example:

   replot.r prof -I k4.test.zip -O k4.replot -SE 0 -MW 9 -H 0.3

   Here I re-create the average profile of H3k4me3 with standard errors 
   turned off. I want to smooth the curve so I use "-MW 9". I also add some 
   shade under the curve for beautifiation by "-H 0.3".

   Another usage with replot.r is that you can explore different gene ranking 
   algorithms in the heatmap. By default, ngs.plot ranks the genes by the toal 
   enrichment of the first profile in each region. However, there many other 
   options you can do. For example, I want to see how hierarchical clustering 
   would group those genes. I can do this:

   replot.r heatmap -I encode1M.k4k27.zip -O k4k27.replot -GO hc -RR 80

   This will perform hierarchical clustering on the genes. I also increase the
   value of reduce ratio (-RR) so that I can get a shorter heatmap because my
   screen is not as large as an IMAX silver screen to contain all genes :)

4. Finally I show you an example of using gene lists in configuration file. In 
   the example folder, there are three example gene lists: "high.Ensembl.txt", 
   "mid.Ensembl.txt" and "low.Ensembl.txt". Each list contains 2,000 Ensembl 
   gene ids that are selected base on gene expression levels derived from RNA-
   seq experiments. The three lists represent genes that are expressed at 
   different levels: high, medium and low. Let's investigate the enrichment of 
   the two histone marks, H3k4me3 and H3k27me3, at the three different 
   expression levels. In order to do that, we create a configuration file:

   hesc.H3k4me3.1M.bam<TAB>high.Ensembl.txt<TAB>"H3k4me3 high"
   hesc.H3k4me3.1M.bam<TAB>mid.Ensembl.txt<TAB>"H3k4me3 mid"
   hesc.H3k4me3.1M.bam<TAB>low.Ensembl.txt<TAB>"H3k4me3 low"
   hesc.H3k27me3.1M.bam<TAB>high.Ensembl.txt<TAB>"H3k27me3 high"
   hesc.H3k27me3.1M.bam<TAB>mid.Ensembl.txt<TAB>"H3k27me3 mid"
   hesc.H3k27me3.1M.bam<TAB>low.Ensembl.txt<TAB>"H3k27me3 low"

   like above and named as "config.genelist.txt". Then we issue a command:

   ngs.plot.r -G hg19 -R tss -C config.genelist.txt -O k4k27.expgrp -D ensembl

   Please note that "-D ensembl" option is used because the gene ids are in 
   Ensembl format. If ngs.plot cannot identify gene ids/symbols with the 
   database, it will emit an error and exit the program.

5. New example: how to use bam pairs? A bam pair is a pair of bam files which 
   are separated by colon. This can be used to represent antibody target vs. 
   DNA or IgG control input. ngs.plot will calculate log2 ratios and plot them 
   using appropriate colors. With the example files, issue a command like this:

   ngs.plot.r -G hg19 -R tss -C hesc.H3k4me3.1M.bam:hesc.Input.500K.bam -O 
   k4vsInp

   Please note that in the heatmap, green colors represent negative values and 
   red colors represent positive values. If you want to create a 
   configuration, do NOT mix bam files with bam pairs because ngs.plot won't 
   be able to plot the two types of curves together.


========Example for RNA-seq data:========

hesc.RNAseq.1M.bam is the alignment file for RNA-seq from ENCODE hESC cell line
. Again, here I only kept 1M reads. Typically we use "genebody" region to draw 
plots for RNA-seq, and we must tell ngs.plot that the input is RNA-seq. This 
will make ngs.plot remove intronic regions to make a correct representation 
for transcripts. To demonstrate, let's issue a command like:

   ngs.plot.r -G hg19 -R genebody -F rnaseq -C hesc.RNAseq.1M.bam -O encode1M.
   rnaseq

Note that an option "-F rnaseq" is used here to turn on the RNA-seq switch.


========Additional notes about gene list files:========

Gene lists can be used to subset genomic coverage. A list can be used directly 
in command line or be supplied in a configuration file. The format is simply 
one gene/transcript symbol/ID per line. Both RefSeq and Ensembl annotations are
accepted. When ngs.plot sees a gene, it actually uses the longest transcript 
of that gene as a representative. You can even mix symbols and IDs and genes 
and transcripts. But be careful that ngs.plot is case-sensitive.


========Extract alternative spliced regions from Cufflinks results:========

Cufflinks must have been run with cuffmerge and cuffcompare to generate a full
set of *.diff files and a GTF file representing the assembled transcripts.
Make sure all files are under the same folder, then issue a command like:

   alt_reg_cufflinks mycuff.gtf res_output

it will generate a bunch of files which contain information such as the
genomic coordinates of alternatively spliced exons, TSS. The BED files can
then be used in ngs.plot to investigate epigenomic changes. A toy example is 
given in "cufflinks_eg" folder. 



For more details of the usage, please visit project website:
http://code.google.com/p/ngsplot/








