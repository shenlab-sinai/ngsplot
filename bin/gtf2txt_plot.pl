#!/usr/bin/perl -w
# Analyze a GTF file downloaded from UCSC or Ensembl, annotate all exons
# and output genebody and exonic coordinates to a TABULAR file for ngs.plot.

use strict;
use MyBioinfo::GTF;

# Program arguments.
if(@ARGV<1){
	print "Usage: $0 database_name input.gtf output.txt\n\n";
	print "database_name    type either refseq or ensembl.\n";
	exit 0;
}

# Database name check.
my $database = $ARGV[0];
unless($database eq 'refseq' or $database eq 'ensembl'){
	print "Type refseq or ensembl for database name. Stop.\n";
	exit 0;
}

# Read gene structures into a hash table.
my $gtf = $ARGV[1];
my %gene_table = read_gtf_tbl($gtf);
anno_gene_table(\%gene_table);

# Output exonic as well as intronic coordinates with annotation.
my $bed = $ARGV[2];
my $beh;
open $beh, ">$bed" or die "Open output BED file error: $!\n";
while(my($gene_id, $gene_struct) = each %gene_table){
	&output_gene($gene_id, $gene_struct, $beh, $database);
}
close $beh;



###### Start Subroutines #######
# Output a whole gene structure.
sub output_gene{
	my($gene_id,$gene_struct,$beh,$database) = @_;
	my $gene_name = $gene_struct->{name};
	while(my($transcript_id,$trans_struct) = each %{$gene_struct->{transcripts}}){
		my $trans_desc;
		if($database eq 'ensembl'){
			$trans_desc = join("\t", $gene_id, $gene_name, $transcript_id);
		}elsif($database eq 'refseq'){
			$trans_desc = join("\t", 'NA', $gene_id, $transcript_id);
		}else{
			die "Unrecognized database name. Exit.\n";
		}
		&output_transcript($gene_struct->{chrom}, $trans_desc, 
			$gene_struct->{strand}, $trans_struct->{exons},$beh);
	}
}

# Output a transcript structure.
# Note: Convert 1-based coordinates in GTF to 0-based in BED!
sub output_transcript{
	my($chrom,$trans_desc,$strand,$exons,$h) = @_;
	# Put the whole transcript as one line in file.
	my $first_exon = $exons->[0];
	my $last_exon = $exons->[$#{$exons}];
	my $trans_start = $first_exon->{'start'};
	my $trans_end = $last_exon->{'end'};
	print $h join("\t", $chrom, $trans_start, $trans_end, $trans_desc, $strand, 'genebody', 'NA'), "\n";
	foreach my $e(@{$exons}){
		# Output exon.
		print $h join("\t", $chrom, $e->{start}, $e->{end}, $trans_desc, $strand, 'exon',  $e->{class}), "\n";
	}
}


# Example lines.
#chr1	pseudogene	exon	3044314	3044814	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "Gm16088"; oId "ENSMUST00000160944"; nearest_ref "ENSMUST00000160944"; class_code "=";
#chr1	snRNA	exon	3092097	3092206	.	+	.	gene_id "XLOC_000002"; transcript_id "TCONS_00000002"; exon_number "1"; gene_name "U6.149"; oId "ENSMUST00000082908"; nearest_ref "ENSMUST00000082908"; class_code "=";
#chr1	processed_transcript	exon	3456668	3456768	.	+	.	gene_id "XLOC_000003"; transcript_id "TCONS_00000003"; exon_number "1"; gene_name "Gm1992"; oId "ENSMUST00000161581"; nearest_ref "ENSMUST00000161581"; class_code "="; tss_id "TSS1";
#chr1	processed_transcript	exon	3503486	3503634	.	+	.	gene_id "XLOC_000003"; transcript_id "TCONS_00000003"; exon_number "2"; gene_name "Gm1992"; oId "ENSMUST00000161581"; nearest_ref "ENSMUST00000161581"; class_code "="; tss_id "TSS1";
#chr1	snRNA	exon	4519098	4519204	.	+	.	gene_id "XLOC_000004"; transcript_id "TCONS_00000004"; exon_number "1"; gene_name "U6.1"; oId "ENSMUST00000082442"; nearest_ref "ENSMUST00000082442"; class_code "=";
#chr1	protein_coding	exon	4797869	4798063	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000005"; exon_number "1"; gene_name "Lypla1"; oId "ENSMUST00000134384"; nearest_ref "ENSMUST00000134384"; class_code "="; tss_id "TSS2";
#chr1	protein_coding	exon	4798536	4798567	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000005"; exon_number "2"; gene_name "Lypla1"; oId "ENSMUST00000134384"; nearest_ref "ENSMUST00000134384"; class_code "="; tss_id "TSS2";
#chr1	protein_coding	exon	4818665	4818730	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000005"; exon_number "3"; gene_name "Lypla1"; oId "ENSMUST00000134384"; nearest_ref "ENSMUST00000134384"; class_code "="; tss_id "TSS2";
#chr1	protein_coding	exon	4820349	4820396	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000005"; exon_number "4"; gene_name "Lypla1"; oId "ENSMUST00000134384"; nearest_ref "ENSMUST00000134384"; class_code "="; tss_id "TSS2";
#chr1	protein_coding	exon	4822392	4822462	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000005"; exon_number "5"; gene_name "Lypla1"; oId "ENSMUST00000134384"; nearest_ref "ENSMUST00000134384"; class_code "="; tss_id "TSS2";
