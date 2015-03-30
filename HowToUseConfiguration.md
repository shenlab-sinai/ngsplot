A config file is composed with one profile per line. Each profile contains fields that are separated by TAB. The whole file may contain 3 to 5 columns. The first 3 columns are required while the other 2 are optional.

The use of gene lists is explained in this wiki: https://code.google.com/p/ngsplot/wiki/HowtoMakeThingsRight OR use -1 for the whole genome. You can also use BED files in this column.

The Bam column may contain a bam file or a bam-pair. Correspondingly, the Fragment Length column may also be a single integer or an integer pair like 150:300.

The color column is used to create custom colors for average profile plot. You should specify R colors here, the same way you specify colors in "plot". If you want to use custom colors without changing the default fragment length, you need to put 150 on the 4th column just to "fill the gap". The color information has to be on the 5th column.

### Table: explanation of a config file ###
| **Bam** | **Gene List** | **Title** | **Fragment Length (Optional)** | **Color (Optional)** |
|:--------|:--------------|:----------|:-------------------------------|:---------------------|
| your.bam | your.txt | "My Sample" | 150 | green |
| your.bam | your.bed | "My Sample2" | 200 | orange |
| your.bam:mine.bam | our.txt | "You vs. Me" | 150:300 | blue |

### Example config files ###
```
# Example 1:
A.bam<TAB>a.txt<TAB>"NameA"
B.bam<TAB>b.txt<TAB>"NameB"
```

```
# Example 2:
A.bam<TAB>a.bed<TAB>"NameA"<TAB>150<TAB>blue2
B.bam<TAB>b.bed<TAB>"NameB"<TAB>300<TAB>green
```

```
# Example 3:
A.bam:B.bam<TAB>-1<TAB>"A vs. B"<TAB>150:300
```