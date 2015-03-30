Bed format is a text-based format that contain 3 required columns: chrom, start, end. ngs.plot can also use the next 3 optional columns: name, score, strand. The name column is interpreted as gene name while the score column is interpreted as transcript id. They are no longer in use from v2.07 because we encourage you to use individual bed files for subset regions. The strand column can be useful since ngs.plot will flip the coverage vector when it sees "-" strands.

A common mistake that people make about bed format is ignoring that bed is: 0-based and right-side open. So the end position is at least 1bp larger than the start position. For example, the following bed line:

```
chr1<TAB>100<TAB>101
```

means the 101st position on chr1 with 1bp in size.