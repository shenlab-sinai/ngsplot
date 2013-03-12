SplitDB <- function(d, s, g){
# Split the original ngs.plot databases into smaller ones based on functions.
# This aims to improve the efficiency in loading genomic coordinates.
# Args:
#   d: database object.
#   s: species name.
#   g: database name, such as refseq or ensembl.
    
	# genebody
    d$genebody -> genome.coord
    genome.coord[order(genome.coord$chrom, genome.coord$start), ] -> genome.coord
    out.f <- sprintf("%s.%s.genebody.RData", s, g)
    save(genome.coord, file=out.f)
   
    require(foreach)
    
	# exon
    foreach(n=names(d$exon)) %do% {
        d$exon[[n]] -> genome.coord
        genome.coord[order(genome.coord$chrom, genome.coord$start), ] -> genome.coord
        out.f <- sprintf("%s.%s.exon.%s.RData", s, g, n)
        save(genome.coord, file=out.f)
    }

	# cgi
    foreach(n=names(d$cgi)) %do% {
        d$cgi[[n]] -> genome.coord
        genome.coord[order(genome.coord$chrom, genome.coord$start), ] -> genome.coord
        out.f <- sprintf("%s.%s.cgi.%s.RData", s, g, n)
        save(genome.coord, file=out.f)
    }
    
	# exonmodel
    d$exonmodel -> exonmodel
    out.f <- sprintf("%s.%s.exonmodel.RData", s, g)
    save(exonmodel, file=out.f)
}



