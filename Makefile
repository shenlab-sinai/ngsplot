BINCUFF=bin/alter2bed.pl bin/combine_diff.pl bin/coordinat.pl bin/difflist2bed.pl bin/get_difflist.pl bin/gtf2tree.pl bin/parse_diff.pl bin/alt_reg_cufflinks
BINNGSP=bin/cov.calc.r bin/ngs.plot.r bin/replot.r
LIB=lib/parse.args.r lib/plotmat.r lib/sep.filename.r lib/smoothplot.r
OTHERS=database example README Changes

ngsplot-dist.tar.gz: $(BINCUFF) $(BINNGSP) $(LIB) $(OTHERS)
	rm -rf ngsplot
	mkdir -p ngsplot/bin ngsplot/lib
	cp -r $(BINCUFF) $(BINNGSP) ngsplot/bin
	cp -r $(LIB) ngsplot/lib
	cp -r $(OTHERS) ngsplot
	tar czvf ngsplot-dist.tar.gz ngsplot
	rm -rf ngsplot
clean:
	rm -rf ngsplot-dist.tar.gz
