# Introduction #
Here is an example workflow which demonstrates the use of ngs.plot to study the enrichment of H3K4me3 at TSS in hesc cell line. The ChIP-seq data is available at GEO: GSM733680.


# Details #
**Home page**
<img src='http://ngsplot.googlecode.com/files/webngsplot1.png' />
```
Figure 1. Home page of ngs.plot Galaxy portal.
```


**Step1: Upload the alignment file to Galaxy (Fig. 2). The BAM file (hesc.H3K4me3.bam) represents the alignment of H3K4me3 ChIP-seq reads.**
<img src='http://ngsplot.googlecode.com/files/webngsplot2.png' />
```
Figure 2. Upload the BAM file - “hesc.H3K4me3.bam” to Galaxy.
```


**Step2: Use ngs.plot to plot the H3K4me3 enrichment at the TSS (Fig. 3). In this page, user can choose the species of the sample. There are 3 model species supported at this moment: human (hg19), mouse (mm9), and rat (rn4). User can also choose parameter value to filter reads which have low mapping quality and expected fragment length. Many graphic options are also provided for user. User can choose from a list of predefined genomic regions, a flanking region size, or specify an
image title. User can also choose to upload a gene list file, a random sampling rate or a smooth method.**
<img src='http://ngsplot.googlecode.com/files/webngsplot3.png' />
```
Figure 3. Use ngs.plot to plot the enrichment at the TSS.
```


**Finally, three files will be produced by ngs.plot: an avgprof file - "AVG file on hesc.H3K4me3.bam", heatmap file - "Heatmap file on hesc.H3K4me3.bam", and a zip file - "Zip file on hesc.H3K4me3.bam"
which can be downloaded for reproducing the results in other programs such as Excel (Figs. 4,5).**
<img src='http://ngsplot.googlecode.com/files/webngsplot4.png' />
```
Figure 4. The enrichment of H3K4me3 at TSS.
```
<img src='http://ngsplot.googlecode.com/files/webngsplot5.png' />
```
Figure 5. The Heat map of H3K4me3 at TSS.
```