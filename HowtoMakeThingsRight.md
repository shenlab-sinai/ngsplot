## Gene List ##
Weird gene names do exist, such as "Olfr1359#2", "Mageb1#2", "CYCA2;1" and "beta'Cop". Because "#" is commonly used as a symbol for comment, the rest part after that on a line will be eaten by a text parser. While "\'" is widely used as quote, and the content between two "\'" will be treated as one column. These will cause the gene names in your list to mismatch with that in the database. To make sure your gene names are properly read, use double quotes. For example, your gene list should look like this:
```
"Gm5512#2"
"Gm5643#2"
"Gm3286#2"
"Mir684-1#4"
"Olfr1359#2"
"beta'Cop"
"CYCA2;1"
```

## Get Ordered Gene List ##
Here is an R script to get the gene names as ordered in the heatmap. Assume you have saved it as ExtractGName.R and your ngs.plot output zip file is: ngsres.zip

```
Usage: Rscript ExtractGName.R ngsres
```

The ordered gene list will be saved in ngsres.gname.txt.

```
fname <- commandArgs(T)

zip.fname <- paste(fname, 'zip', sep='.')
heatmap.dat <- file.path(fname, 'heatmap.RData')
load(unz(zip.fname, heatmap.dat))

gene.list <- strsplit(go.list[[1]], ':')
gname.list <- sapply(gene.list, function(x) x[1])
write.table(gname.list, file=paste(fname, 'gname.txt', sep='.'),
            col.names=F, row.names=F)
```