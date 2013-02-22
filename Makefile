BINCUFF=bin/alter2bed.pl bin/combine_diff.pl bin/coordinat.pl bin/difflist2bed.pl bin/get_difflist.pl bin/gtf2tree.pl bin/parse_diff.pl bin/alt_reg_cufflinks
BINNGSP=bin/cov.calc.r bin/ngs.plot.r bin/replot.r
LIB=lib/parse.args.r lib/plotlib.r lib/sep.filename.r lib/coverage.r lib/genedb.r
ALLR=${BINNGSP} ${LIB}
OTHERS=database example README Changes
DISTFOLDER=ngsplot

ngsplot-dist.tar.gz: $(BINCUFF) $(BINNGSP) $(LIB) $(OTHERS)
	rm -rf ${DISTFOLDER}
	mkdir -p ${DISTFOLDER}/bin ${DISTFOLDER}/lib
	cp -r $(BINCUFF) $(BINNGSP) ${DISTFOLDER}/bin
	cp -r $(LIB) ${DISTFOLDER}/lib
	cp -r $(OTHERS) ${DISTFOLDER}
	# Remove comments and blank lines from *.r files.
	for r in ${ALLR}; do \
		sed '2,+10000 s/^\s*#.*//' ${DISTFOLDER}/$$r|sed '/^$$/ c\' > tmp; \
		mv tmp ${DISTFOLDER}/$$r; \
	done
	tar czvf $@ ${DISTFOLDER}
	rm -rf ${DISTFOLDER}

clean:
	rm -rf ngsplot-dist.tar.gz
