# Arguments specific to ngs.plot.r #
## Mandatory arguments ##
| **Argument** | **Explanation** | **Accepted value and notes** |
|:-------------|:----------------|:-----------------------------|
| **-G** | Genome name | currently supported: see wiki page ("ID" column) https://code.google.com/p/ngsplot/wiki/SupportedGenomes |
| **-R** | Genomic regions to plot | tss, tes, genebody, exon, cgi, enhancer, dhs or bed(custom regions) |
| **-C** | Bam file or a configuration file for multiple plot | Configuration file must end with ".txt" (See config.example.txt). A new function allows you to specify a bam pair like: A.bam:B.bam. A bam-pair can be used to represent antibody target vs. DNA input or IgG control. ngs.plot will calculate log2 ratios and plot them using appropriate colors. |
| **-O** | Name of output | Multiple files will be generated with appropriate suffixes added, including two PDFs of average profile and heatmap and a zip file of result data. |

## Parameters that can be provided here or in configuration file for multiplot ##
| **Argument** | **Explanation** | **Accepted value and notes** |
|:-------------|:----------------|:-----------------------------|
| **-E** | Gene list to subset regions OR bed file for custom region | If no list is given, the whole genome will be plotted. A text file should be provided with each gene ID or symbol per line. Both gene symbols and Ensembl IDs are accepted. You can also provide transcript IDs and even mix them with gene symbols/IDs. When use -R bed, you can also supply your own bed file here. |
| **-T** | Image title | The title will be shown in the figure legend. Default=Noname.|

## Coverage-generation parameters ##
| **Argument** | **Explanation** | **Accepted value and notes** |
|:-------------|:----------------|:-----------------------------|
| **-F** | Further information provided to select database table or plottype | This is a string of description separated by comma. E.g. protein\_coding,K562,rnaseq (order of descriptors does not matter) means coding genes in K562 cell line drawn in rnaseq mode. See Wiki page: https://code.google.com/p/ngsplot/wiki/UseFurtherInfo |
| **-D** | Gene database | ensembl (default), refseq. |
| **-L** | Flanking region size | Size in bps. By default, when -R=tss, tes, genebody, -L=2000; when -R=exon, cgi, -L=500; when -R=`*`.bed, -L=1000. Will override flanking region factor. |
| **-N** | Flanking region factor | When chosen, the flanking region size equals to the interval size times flanking factor. It can be any numeric value. This allows the flanking region size to change dynamically and may create plots that look more "natural". |
| **-RB** | Robust statistics fraction | The fraction of extreme values equal to this amount will be trimmed on both ends. By default, this is set to 0 (0%). Use 0.05 to remove 5% of extreme values.|
| **-S** | Random sampling rate | Must be (0, 1]. This will randomly sample a portion of the whole genome or the gene lists. This option can be VERY useful if one just wants to get an overview in shorter time. |
| **-P** | #CPUs to be used | Set 0 to use all CPUs that are detected on your machine. |
| **-AL** | Algorithm used to normalize coverage vectors: | Coverage vectors can be of various lengths, such as those from genebodies. They must be normalized to equal length to be averaged and plotted. |
|  | **spline(default)** | Spline fit first and then values are taken at equal intervals. |
|  | **bin** | The whole vector is split into a fixed number of equal sized bins and the average for each bin is calculated. |
| **-CS** | Chunk size for loading genes in batch. | This parameter controls the behavior of coverage calculation. A smaller value implies lower memory footprint but may increase processing time. I suggest to leave it as is unless you really understand what you are doing! |
| **-MQ** | Mapping quality cutoff to filter reads | Default=20. This is the Phred-scale mapping quality score. A score of 20 means an error rate of 1%. |
| **-FL** | Fragment(insert) length used to calculate physical coverage | Default=100. ngs.plot calculates physical instead of read coverage. This will produce figures that contain more accurate representation of ChIP enrichment. You should set this value equal to the average fragment length in your sequencing library. |
| **-SS** | Strand-specific coverage calculation | Choose: both(default), same, opposite. By default, the short reads from both strands will be counted. However, you may specify if you want only the reads that are on the same or the opposite strand as the transcript. |
| **-IN** | Boolean tag for large interval | 0 or 1. By default, exon and cgi are small interval; genebody and `*`.bed are large interval. The X-axis of the plot is separated into 5 equal sized parts. If choose large interval, the middle 3 parts are for interval region and the 1 part on each side is for flanking region. If choose small interval, the middle 1 part is for interval region and the 2 parts on each side is for flanking region. |
| **-FI** | Boolean tag for forbidding image output | 0 or 1. Set 1 to turn off image output. If you are running ngs.plot on a remote server which does not support graphic output, this can be useful. Transfer the text output to your local machine and use replot.r to generate images. |


---

# Arguments common to both ngs.plot.r and replot.r #
## Miscellaneous parameters ##
| **Argument** | **Explanation** | **Accepted value and notes** |
|:-------------|:----------------|:-----------------------------|
| **-FS** | Font size | Default is 12 pt. |

## Avg. profile parameters ##
| **Argument** | **Explanation** | **Accepted value and notes** |
|:-------------|:----------------|:-----------------------------|
| **-WD** | Image width | Default is 8 in. |
| **-HG** | Image height | Default is 7 in. |
| **-SE** | Boolean tag for plotting standard errors | 0 or 1. By default, standard errors will be rendered as shaded area around each curve. |
| **-MW** | Moving window width to smooth avg. profiles | Default=1(no smoothing). The unit of the window size is a data point. In ngs.plot, the number of data points of the X-axis is 100. A moving window width of 10 means averaging over 10% of the entire X-axis. |
| **-H** | Opacity of shaded area | Suggested value: [0, 0.5]. This will add semi-transparent shade under each curve. |
| **-YAS** | Y-axis scale | auto(default) or min\_val,max\_val(custom scale). |
| **-LEG** | Draw legend? | Boolean to control the display of legend: 1(default) or 0. |
| **-BOX** | Draw box around plot? | Boolean to control the display of box: 1(default) or 0. |
| **-VLN** | Draw vertical lines? | Boolean to control the display of vertical lines(e.g., TSS and TES): 1(default) or 0. |
| **-XYL** | Draw X- and Y-axis labels? | Boolean to control the display of X- and Y-axis labels: 1(default) or 0. |
| **-LWD** | Line width | Default is 3 pt. |


## Heatmap parameters ##
| **Argument** | **Explanation** | **Accepted value and notes** |
|:-------------|:----------------|:-----------------------------|
| **-GO** | Gene order algorithm | The algorithm is used to rank genes in heatmaps. Available options: |
|  | **total(default)** | Overall enrichment of the 1st profile |
|  | **hc** | Hierarchical clustering |
|  | **max** | Peak value of the 1st profile. This option makes more sense if the epigenomic mark tends to generate sharper peak. |
|  | **prod** | Product of all profiles on the same region. |
|  | **diff** | Difference between the 1st and the 2nd profiles. |
|  | **km** | K-means clustering. The default number of clusters is 5. |
|  | **none** | No ranking algorithm applied. Use order provided in the gene list. This can be used to your advantage. For example, you can rank genes by their expression levels and give the list to ngs.plot and see how the epigenomes change with expression. |
| **-LOW** | Low count cutoff in rank-based normalization. | This is raw read count. The default is 10. |
| **-KNC** | K-means number of clusters | The default is 5. |
| **-MIT** | Maximum number of iterations for K-means | The default is 20. |
| **-NRS** | Number of random starts in K-means. | K-means is prone to local optima. Restarting it repeatedly may help to find a better solution. The default is 30. |
| **-RR** | Reduce ratio | The parameter controls the heatmap height. The smaller the value, the taller the heatmap. The default is 30. |
| **-SC** | Color scale for heatmaps | This sets a range for mapping data values to colors. Any values belong this range will be mapped to the same colors as the range ends. |
|  | **local(default)** | Each individual heatmap has its own color scale. |
|  | **region** | All heatmaps that belong to the same region have the same color scale. |
|  | **global** | All heatmaps in the current plot have the same color scale. |
|  | **min\_val,max\_val** | Custom color scale using a pair of numerics. E.g., 0,5 means the smallest value is 0 and the largest value is 5. |
| **-FC** | Flooding fraction | Default=0.02(2%). That means, the minimum value is truncated at 2% and the maximum value is truncated at 98%. A higher fraction results in plots that have higher brightness but are less dynamic. |
| **-CO** | Color for heatmap | For bam-pair, use color-tri(neg\_color:`[`neu\_color`]`:pos\_color). Hint: must use R colors, such as darkgreen, yellow and blue2. The neutral color is optional(default=black). |
| **-CD** | Color distribution for heatmap | The default is 0.6. Must be a positive number. Hint: lower values give more widely spaced colors at the negative end. In other words, they shift the neutral color to positive values. If set to 1, the neutral color represents 0(i.e. no bias). If set to >1, the neutral color represents a negative value. |
