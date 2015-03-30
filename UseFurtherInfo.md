# Introduction #

ngs.plot now contains a database of thousands of tables across multiple species, gene types and regions. This makes it difficult for a user to find the genomic regions for plot. To avoid cluttering the command line interface, we design a mechanism for users to select tables based on a single string of descriptors that are separated by comma. This keeps the command line concise while remains the flexibility.


# Usage #

Briefly, the "-F" option should be followed by a string of descriptors like this:

```
-F [gene_type][,sub_region][,cell_line or tissue][,exon_type][,rnaseq or chipseq]
```

Please be advised that the option to switch between chip-seq and rna-seq plots can also be provided here. You don't have to give all the descriptors. One descriptor or three descriptors work just fine. When a certain variable is not specified, ngs.plot chooses the default values. For example,

```
-G hg19 -R enhancer -F K562
```

will select the enhancers in the K562 cell line for protein coding genes, although "protein\_coding" is not specified here. The order of the descriptors does not matter, that means

```
-F protein_coding,K562
```

and

```
-F K562,protein_coding
```

are the same. Do NOT feel nervous for making typos or inputting wrong combinations. ngs.plot will always choose a default table if some descriptions do not match the database. The actual database table being used will be given on the command console for your reference.

# Descriptors #

Here are the available values you can use in ngs.plot:

**Gene Type**
|protein\_coding(default) |pseudogene |lincRNA |miRNA |misc(everything else)|
|:------------------------|:----------|:-------|:-----|:--------------------|

**Exon Type**
|canonical(default) |altAcceptor |altBoth |altDonor |polyA |promoter |variant|
|:------------------|:-----------|:-------|:--------|:-----|:--------|:------|

**CGI & DHS classification**
|ProximalPromoter(default) |Promoter1k |Promoter3k |Genebody |Genedesert |OtherIntergenic |Pericentromere |Subtelomere|
|:-------------------------|:----------|:----------|:--------|:----------|:---------------|:--------------|:----------|

**Human cell line**

dhs

|H1hesc(default) |A549 |Gm12878 |Helas3 |Hepg2 |Hmec |Hsmm |Hsmmtube |Huvec |K562 |Lncap |Mcf7 |Nhek |Th1|
|:---------------|:----|:-------|:------|:-----|:----|:----|:--------|:-----|:----|:-----|:----|:----|:--|

enhancer

|H1hesc(default) |Gm12878 |Hepg2 |Hmec |Hsmm |Huvec |K562 |Nhek |Nhlf|
|:---------------|:-------|:-----|:----|:----|:-----|:----|:----|:---|

**Mouse tissue(enhancer)**
|mESC(default) |boneMarrow |cerebellum |cortex |heart |intestine |kidney |liver |lung |MEF |olfactoryBulb |placenta |spleen |testes |thymus|
|:-------------|:----------|:----------|:------|:-----|:---------|:------|:-----|:----|:---|:-------------|:--------|:------|:------|:-----|