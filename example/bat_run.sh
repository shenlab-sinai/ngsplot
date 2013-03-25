ngs.plot.r -G hg19 -R tss -C hesc.H3k4me3.1M.bam -O k4.test
ngs.plot.r -G hg19 -R tss -C config.example.txt -O encode1M.k4k27
replot.r prof -I k4.test.zip -O k4.replot -SE 0 -MW 9 -H 0.3
replot.r heatmap -I encode1M.k4k27.zip -O k4k27.replot -GO hc -RR 80
ngs.plot.r -G hg19 -R genebody -F rnaseq -C hesc.RNAseq.1M.bam -O encode1M.rnaseq
ngs.plot.r -G hg19 -R tss -C hesc.H3k4me3.1M.bam:hesc.Input.500K.bam -O k4vsInp
