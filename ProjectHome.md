# INTRODUCTION #
**ngs.plot is a program that allows you to easily visualize your next-generation sequencing (NGS) samples at functional genomic regions.**

DNA sequencing is at the core of genomics. The NGS technology has been tremendously improved in the past few years. It can now determine more than a billion DNA sequences within a week, generating terabytes of data. Applications include but are not limited to: 1. ChIP-seq which profiles genome-wide protein-DNA interactions; 2. RNA-seq which measures the gene expression levels. It is very helpful to look at the enrichment of those sequences at various functional regions. Although a genome browser (such as the UCSC genome browser) allows a researcher to visualize these data, it limits the view to a slice of the genome. While the genome is like a huge collection of functional elements that can be classified into different categories. Each category of elements may perform distinct functions and they might further contain modules.

The signature advantage of ngs.plot is that it collects a large database of functional elements for many genomes. A user can ask for a functionally important region to be displayed in one command. It handles large sequencing data efficiently and has only modest memory requirement. For example, ngs.plot was used to draw a plot for all the genes on the mouse genome from 71GB of ChIP-seq data in 25 min, with a memory footprint of 2.7GB using 4 x 2.4GHz CPU cores. ngs.plot is also easy to use. A user only needs to create a very small text file called configuration, telling the program which samples to look at and how they should be combined with different regions, and then run the program with one command. A web-based version (integrated into Galaxy) is also available for the ones who are allergic to terminals.

## Program Download Location ##
Since Google canceled the download function on Google code, we have to move all of our future program download files to this Google drive folder:

https://drive.google.com/folderview?id=0B1PVLadG_dCKN1liNFY0MVM1Ulk&usp=sharing

**The last release on the download page is v2.08**. For more recent releases, please visit the above link.

## Supported Genomes ##
ngs.plot has an approach to install genomes on demand. It can support for any genome. All you need to do is to download an archive file and install it by yourself. The genome files can be found in this Google drive folder:

https://drive.google.com/folderview?id=0B1PVLadG_dCKNEsybkh5TE9XZ1E&usp=sharing

If you cannot find yours, please request your genome of interest in this post:

https://groups.google.com/forum/#!topic/ngsplot-discuss/GyoacAzV7jM

A list of the available genomes is also documented in this Wiki-page: https://code.google.com/p/ngsplot/wiki/SupportedGenomes

A brief list is here: human (hg18, hg19), chimpanzee (panTro4), rhesus macaque (rheMac2), mouse (mm9, mm10), rat (rn4, rn5), cow (bosTau6), chicken (galGal4), zebrafish (Zv9), drosophila (dm3), Caenorhabditis elegans (ce6, ceX), Saccharomyces cerevisiae (sacCer2, sacCer3), Schizosaccharomyces pombe (Asm294), Arabidopsis thaliana (TAIR10), Zea mays (AGPv3), rice (IRGSP-1.0).

## Extension Annotation Package: Enhancers and DHSs ##
There are now two ngs.plot databases for the reference genome hg19 and mm9. The basic package is a light version of the genome, which only contains the regular stuffs like genebody, CGI and exon. While extension package is an extension of the genome, which contains enhancers or dhs.

By default, the main program installation comes with basic package only. However, if you want to plot enhancers and dhs for the two genomes, you can download the package files from the shared Google drive folder and install them using ngsplotdb.py.

## A note about the compatibility between genome files and program versions ##
The 3.0+ genome files work with ngs.plot v2.41+; The 1-2 series of genome files work with ngs.plot v2.08 and below.

# RECENT CHANGES #
Recent changes, bug fixes and feature additions will be announced through this Google discussion group: https://groups.google.com/forum/?fromgroups#!forum/ngsplot-discuss

If you are interested, please sign up to receive updates through E-mails.

# BLEEDING EDGE #
ngs.plot development has been migrated to Google code's git server. You can checkout a clone at the "Source" tab. This method allows us to respond more quickly to issues than new package releases. Users should always assume the newest version is at the develop branch of the git repo. If you want to become a contributor to ngs.plot, please let me know!

# INSTALLATION #

You need R version >= 2.15.0 and Python 2.7 to be able to use ngs.plot. Please also update related R and Bioconductor packages to the corresponding version.

1. Download the ngs.plot package to a desired folder, such as ~/software, and extract it:
```
cd ~/software
tar xzvf ngsplot-XXX.tar.gz
```
The program will be extracted into a folder called "ngsplot".

2. Add ngsplot executables to your PATH. Under bash, add line like this:
```
export PATH=~/software/ngsplot/bin:$PATH
```
to your ~/.bash\_profile

3. Set environment variable NGSPLOT like this in your ~/.bash\_profile:
```
export NGSPLOT=~/software/ngsplot
```
Then in the terminal, execute:
```
source ~/.bash_profile
```

A trick to avoid setting environment variable is:
```
NGSPLOT=/your/path/to/ngsplot bash -c 'ngs.plot.r XXX'
```

4. Install some required libraries in R:
```
install.packages("doMC", dep=T)
install.packages("caTools", dep=T)
install.packages("utils", dep=T)
```

For R 3.0+, "utils" is no longer needed. For R <3.0, you probably already have caTools and utils
installed but it does not hurt to check.

Then execute in R:
```
source("http://bioconductor.org/biocLite.R")
biocLite( "BSgenome" )
biocLite( "Rsamtools" )
biocLite( "ShortRead" )
```

5. Install ngsplot package in Galaxy:
```
1) Go to $NGSPLOT/galaxy.

2) Copy "ngsplot", "multiplot" and "replot" folders to <your galaxy_install_dir>/tools/ directory. 

3) Edit  <your galaxy_install_dir>/tool_conf.xml and add the following lines in the <toolbox> section:
   <section name="NGSPLOT Tools" id="ngs_plot">
    <tool file="ngsplot/ngs.plot.xml"/>
    <tool file="replot/replot.xml"/>
    <tool file="multiplot/multiplot.xml"/>
   </section>

4) Restart Galaxy.
```

Here is an example workflow which demonstrates the use of ngs.plot in Galaxy:
http://code.google.com/p/ngsplot/wiki/webngsplot


# USAGE #

A wiki-page has been created for detailed explanation of each argument: http://code.google.com/p/ngsplot/wiki/ProgramArguments101

When you type one of the commands without specifying any argument, the program will print out a brief usage. The following are basically copy & pasted from command console and are provided for quick reference.

<h3>0. NEW to v2.00+</h3> Manipulate annotation database and the use of option "-F" <br>

The ngsplotdb.py script should be easy to use. Here are a few examples:<br>
<br>
<pre><code>ngsplotdb.py list  # List installed genomes.<br>
ngsplotdb.py install ngsplotdb_hg19_71_2.0.tar.gz  # Install reference genome from a package file.<br>
ngsplotdb.py remove hg19  # Remove installed genome.<br>
ngsplotdb.py remove --ftr enhancer hg19  # Remove enhancer installation from hg19.<br>
</code></pre>

The "-F" option is a string of descriptors to refine the annotation to use. See this Wiki page for detailed explanation: <a href='https://code.google.com/p/ngsplot/wiki/UseFurtherInfo'>https://code.google.com/p/ngsplot/wiki/UseFurtherInfo</a>. Here are a few examples:<br>
<br>
<pre><code>-F K562  # Select cell line.<br>
-F K562,lincRNA  # Select cell line and gene type.<br>
-F lincRNA,K562  # Same as above(order does not matter).<br>
-F Promoter3k,H1hesc,protein_coding  # Select region, cell line and gene type(apply to DHS only).<br>
</code></pre>

<h3>I. ngs.plot.r</h3> Use ngs.plot.r to choose a genomic region of interest and create enrichment plots of any ChIP-seq or RNA-seq samples. <br>

<pre><code>Type "ngs.plot.r 2&gt;&amp;1|less" at console will give you a brief usage summary for online <br>
reference. Here is just a truncated output:<br>
<br>
Usage: ngs.plot.r -G genome -R region -C [cov|config]file<br>
                  -O name [Options]<br>
<br>
## Mandatory parameters:<br>
  -G   Genome name. Use ngsplotdb.py list to show available genomes.<br>
  -R   Genomic regions to plot: tss, tes, genebody, exon, cgi, enhancer, dhs or bed<br>
  -C   Indexed bam file or a configuration file for multiplot<br>
  -O   Name for output: multiple files will be generated<br>
</code></pre>

<h3>II. replot.r</h3> Use replot.r to re-create plots with generated data. There are a couple of options for you to fine-tune the figures.<br>

<pre><code>Usage: replot.r command -I input.zip -O name<br>
<br>
  command: prof OR heatmap<br>
<br>
## Mandatory parameters:<br>
    -I  Result zip file created by ngs.plot<br>
    -O  Output name<br>
</code></pre>

<h3>III. plotCorrGram.r</h3> Create a corrgram from ngs.plot output files.<br>

<pre><code>Usage: plotCorrGram.r -I ngsplot_output.zip -O output_name [Options]<br>
<br>
## Mandatory parameters:<br>
  -I   Result zip file created by ngs.plot.<br>
  -O   Output name<br>
## Optional parameters:<br>
  -M   Method used to calculate row stat.<br>
       mean(default): mean of each row.<br>
       max: max of each row.<br>
       window: mean on center region.<br>
  -P   Options for -M method.<br>
       mean: [0,0.5) - trim value for robust estimation, default is 0.<br>
       window: [0,0.5),(0.5,1] - window borders, default:0.33,0.66.<br>
  -D   Options for distance calculation in hierarchical cluster.<br>
       This must be one of 'euclidean'(default), 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.<br>
  -H   Options for agglomeration method in hierarchical cluster.<br>
       This must be one of 'ward'(default), 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.<br>
</code></pre>

<h1>EXAMPLES</h1>


1. ngs.plot.r needs an indexed bam file or a configuration file as an input to plot short read coverage across the genomic regions of interest. ngs.plot.r will generate multiple files including avgprof file, heatmap file, and a zip file for replotting.<br>
<br>
Command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R tss -C hesc.H3k4me3.rmdup.sort.bam -O hesc.H3k4me3.tss -T H3K4me3 -L 3000 -FL 300<br>
</code></pre>

Data from: Ernst, J., Kheradpour, P., Mikkelsen, T.S., Shoresh, N., Ward, L.D., Epstein, C.B., Zhang, X., Wang, L., Issner, R., Coyne, M., et al. (2011). Mapping and analysis of chromatin state dynamics in nine human cell types. Nature 473, 43-49.<br>
<br>
The avgprof and heatmap plotted by ngs.plot are like this:<br>
<img src='http://ngsplot.googlecode.com/files/hesc.H3k4me3.tss.all.png' />

ngs.plot can also accept bam pairs for plot. A bam pair is a pair of bam files separated by colon, such as ChIP vs. Input. Using H3K4me3 as an example, you can issue a command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R tss -C hesc.H3k4me3.rmdup.sort.bam:hesc.Input.rmdup.sort.bam -O hesc.H3k4me3vsInp.tss -T H3K4me3 -L 3000 -FL 300<br>
</code></pre>

The avgprof and heatmap plotted by ngs.plot are like this:<br>
<img src='https://ngsplot.googlecode.com/files/k4_bampair.png' />

2. ngs.plot for multiplot<br>
<br>
If you want to draw a multiplot, you need to configure a txt file for ngs.plot.<br>
<h5><code>(1) H3K4me3</code></h5>
The configure file "config.hesc.k4.txt" is like this:<br>
<pre><code># If you want to specify the gene list as "genome", use "-1".<br>
# Use TAB to separate the three columns: coverage file&lt;TAB&gt;gene list&lt;TAB&gt;title<br>
# "title" will be shown in the figure's legend.<br>
hesc.H3k4me3.rmdup.sort.bam     high_expressed_genes.txt         "High"<br>
hesc.H3k4me3.rmdup.sort.bam     medium_expressed_genes.txt       "Med"<br>
hesc.H3k4me3.rmdup.sort.bam     low_expressed_genes.txt          "Low"<br>
</code></pre>

Command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R genebody -C config.hesc.k4.txt -O hesc.k4.genebody -D ensembl -FL 300<br>
</code></pre>

Data from: ENCODE Project Consortium, et al. (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74.<br>
<br>
The avgprof and heatmap plotted by ngs.plot like this:<br>
<img src='http://ngsplot.googlecode.com/files/hesc.k4.genebody.all.png' />

<h5><code>(2) H3K36me3</code></h5>
The configure file "config.hesc.k36.txt" is like this:<br>
<pre><code># If you want to specify the gene list as "genome", use "-1".<br>
# Use TAB to separate the three columns: coverage file&lt;TAB&gt;gene list&lt;TAB&gt;title<br>
# "title" will be shown in the figure's legend.<br>
hesc.H3k36me3.rmdup.sort.bam     high_expressed_genes.txt         "High"<br>
hesc.H3k36me3.rmdup.sort.bam     medium_expressed_genes.txt       "Med"<br>
hesc.H3k36me3.rmdup.sort.bam     low_expressed_genes.txt          "Low"<br>
</code></pre>

Command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R genebody -C config.hesc.k36.txt -O hesc.k36.genebody -D ensembl -FL 300<br>
</code></pre>

Data from: ENCODE Project Consortium, et al. (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74.<br>
<br>
The avgprof and heatmap plotted by ngs.plot like this:<br>
<img src='http://ngsplot.googlecode.com/files/hesc.k36.genebody.all.png' />

<h5><code>(3) H3K9me3</code></h5>
The configure file "config.hesc.k9.txt" is like this:<br>
<pre><code># If you want to specify the gene list as "genome", use "-1".<br>
# Use TAB to separate the three columns: coverage file&lt;TAB&gt;gene list&lt;TAB&gt;title<br>
# "title" will be shown in the figure's legend.<br>
hesc.H3k9me3.rmdup.sort.bam     high_expressed_genes.txt         "High"<br>
hesc.H3k9me3.rmdup.sort.bam     medium_expressed_genes.txt       "Med"<br>
hesc.H3k9me3.rmdup.sort.bam     low_expressed_genes.txt          "Low"<br>
</code></pre>

Command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R genebody -C config.hesc.k9.txt -O hesc.k9.genebody -D ensembl -FL 300<br>
</code></pre>

Data from: ENCODE Project Consortium, et al. (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74.<br>
<br>
The avgprof and heatmap plotted by ngs.plot like this:<br>
<img src='http://ngsplot.googlecode.com/files/hesc.k9.genebody.all.png' />

<h5><code>(4) H3K27me3</code></h5>
The configure file "config.hesc.k27.txt" is like this:<br>
<pre><code># If you want to specify the gene list as "genome", use "-1".<br>
# Use TAB to separate the three columns: coverage file&lt;TAB&gt;gene list&lt;TAB&gt;title<br>
# "title" will be shown in the figure's legend.<br>
hesc.H3k27me3.rmdup.sort.bam     high_expressed_genes.txt         "High"<br>
hesc.H3k27me3.rmdup.sort.bam     medium_expressed_genes.txt       "Med"<br>
hesc.H3k27me3.rmdup.sort.bam     low_expressed_genes.txt          "Low"<br>
</code></pre>

Command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R genebody -C config.hesc.k27.txt -O hesc.k27.genebody -D ensembl -FL 300<br>
</code></pre>

Data from: ENCODE Project Consortium, et al. (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74.<br>
<br>
The avgprof and heatmap plotted by ngs.plot like this:<br>
<img src='http://ngsplot.googlecode.com/files/hesc.k27.genebody.all.png' />


3. ngs.plot.r can be used to analyze RNA-seq data. Here We used an in-house RNA-seq dataset (unpublished) from human post-mortem brain tissue of schizophrenia patients as an example. One sample with acceptable RNA quality (RIN=7.8) and another sample with degraded RNA quality (RIN=3) are chosen.<br>
<br>
The configure file "config.RIN_number.txt" is like this:<br>
<pre><code># If you want to specify the gene list as "genome", use "-1".<br>
# Use TAB to separate the three columns: coverage file&lt;TAB&gt;gene list&lt;TAB&gt;title<br>
# "title" will be shown in the figure's legend.<br>
Individual1_3.bam                    -1                "Individual1_3"<br>
Individual2_7.8.bam                  -1                "Individual2_7.8"<br>
</code></pre>

Command like this:<br>
<pre><code>ngs.plot.r -G hg19 -R genebody -C config.RIN_number.txt -O RIN_number -F rnaseq<br>
</code></pre>

The avgprof and heatmap plotted by ngs.plot like this:<br>
<img src='http://ngsplot.googlecode.com/files/RIN_number.all.png' />

The plot above shows that the sample with lower RIN number is significantly biased in short read coverage towards the 3’ end.<br>
<br>
<br>
<h2>COMMERCIAL USE</h2>

ngs.plot is free for use by academic users. If you want to use it in commercial settings, please contact Lisa Placanica at: lisa.placanica@mssm.edu<br>
<br>
<br>
<h2>HOW TO CITE</h2>

<b>Shen, L.<code>*</code>, Shao, N., Liu, X. and Nestler, E. (2014) ngs.plot: Quick mining and visualization of next-generation sequencing data by integrating genomic databases, BMC Genomics, 15, 284.</b>


<h2>CONTACT</h2>

ngs.plot is developed by Drs. Li Shen, Ningyi Shao and Xiaochuan Liu at the Icahn School of Medicine at Mount Sinai. If you have technical questions about ngs.plot, please use the discussion forum. For collaborations or any other matters, use: li.shen <b>AT</b> mssm.edu.<br>
