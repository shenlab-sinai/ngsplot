BINCUFF=bin/alter2bed.pl bin/combine_diff.pl bin/coordinat.pl bin/difflist2bed.pl bin/get_difflist.pl bin/gtf2tree.pl bin/parse_diff.pl bin/alt_reg_cufflinks
BINNGSP=bin/cov.calc.r bin/ngs.plot.r bin/replot.r
LIB=lib/parse.args.r lib/plotmat.r lib/sep.filename.r lib/smoothplot.r
OTHERS=database example README Changes

ngsplot.tar.gz: 
	tar czvf ngsplot.tar.gz $(BINCUFF) $(BINNGSP) $(LIB) $(OTHERS)
clean:
	rm -rf ngsplot.tar.gz
