# PERL scripts to extract regions from Cufflinks.
BINCUFF=bin/alter2bed.pl bin/combine_diff.pl bin/coordinat.pl \
		bin/difflist2bed.pl bin/get_difflist.pl bin/gtf2tree.pl \
		bin/parse_diff.pl bin/alt_reg_cufflinks

# ngs.plot R executable R scripts.
BINNGSP=bin/ngs.plot.r bin/replot.r

# ngs.plot R libs.
LIB=lib/parse.args.r lib/plotlib.r lib/coverage.r lib/genedb.r

# All ngs.plot R files.
ALLR=${BINNGSP} ${LIB}

# ngs.plot database, example, readme, etc.
OTHERS=database example README Changes galaxy

# ngs.plot distribution folder
DISTFOLDER=ngsplot

# example bam files.
BAM=example.bam

all: program bam

program: ngsplot-dist.tar.gz
bam: example.bam.tar.gz

ngsplot-dist.tar.gz: $(BINCUFF) $(BINNGSP) $(LIB) $(OTHERS)
	rm -rf ${DISTFOLDER}
	mkdir -p ${DISTFOLDER}/bin ${DISTFOLDER}/lib
	cp -rL $(BINCUFF) $(BINNGSP) ${DISTFOLDER}/bin
	cp -rL $(LIB) ${DISTFOLDER}/lib
	cp -rL $(OTHERS) ${DISTFOLDER}
	find ${DISTFOLDER}/ -name '.svn'|xargs -I% rm -r %
	# Remove comments and blank lines from *.r files.
	# for r in ${ALLR}; do \
	# 	sed '2,+10000 s/^\s*#.*//' ${DISTFOLDER}/$$r|sed '/^$$/ c\' > tmp; \
	# 	mv tmp ${DISTFOLDER}/$$r; \
	# done
	tar czvf $@ ${DISTFOLDER}
	rm -rf ${DISTFOLDER}

example.bam.tar.gz: $(BAM)
	tar czvf $@ $(BAM)

clean:
	rm -rf ngsplot-dist.tar.gz example.bam.tar.gz
