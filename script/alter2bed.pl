#!/usr/bin/perl -w
##
#########################################
##
## File Name:
##   alter2bed.pl
##
## Description:
##   Parse the bed file and get the bed annotation for alternative splicing list.
##
## Usage:
##   $0 <cufflinks_bed_file> <splicing.diff>  <isoform_exp.diff> <alternative_splicing_type> <alternative_splicing.bed> 
##
## Author:
##   Dr. Xiaochuan Liu (xiaochuan.liu@mssm.edu)
##
## Date:
##   Mon Dec 12 09:00:00 EST 2011
##
#########################################
#
#

use strict;
my $usage = "$0 <cufflinks_bed_file> <splicing.diff>  <isoform_exp.diff> <alternative_splicing_type>  <alternative_splicing.bed> \n\nalternative_splicing_type: \n\t1) variant\n\t2) altBoth\n\t3) altDonor\n\t4) altAcceptor\n";

if(@ARGV<5){
        print "Usage: $0  <cufflinks_bed_file> <splicing.diff>  <isoform_exp.diff> <alternative_splicing_type>  <alternative_splicing.bed>\n\nalternative_splicing_type: \n\t1) variant\n\t2) altBoth\n\t3) altDonor\n\t4) altAcceptor\n";
        exit 0;
}


my $cufflinks_bed_file = $ARGV[0];
my $splicing_file = $ARGV[1];
my $isoform_file = $ARGV[2];
my $alternative_splicing_type = $ARGV[3];
my $output_bedfile=$ARGV[4];
$alternative_splicing_type = "exon:".$alternative_splicing_type;

open(IN1, $cufflinks_bed_file) || die $usage ;
open(IN2, $splicing_file) || die $usage ;
open(IN3, $isoform_file) || die $usage ;

my $output_bed;
my %hash_gene_tss = ();
my %hash_gene_trans = ();

open $output_bed, ">$output_bedfile" or die "Open output alternative splicing bed file error: $!\n";

my @cufflinks_array=<IN1>;
my @splicing_array=<IN2>;
my @isoform_array=<IN3>;


for (my $j=0; $j<scalar(@splicing_array); $j++ )
{
     my $splicing_array_string=$splicing_array[$j];
     $splicing_array_string=~s/^\s*|\s*$//g;    #like the String.trim() method
     my @string_splicing_array = split('\t', $splicing_array_string);
     my $tss_id = $string_splicing_array[0];
     my $gene_id = $string_splicing_array[1];
     my $gene_tss = $gene_id.":".$tss_id;
     $hash_gene_tss{$gene_tss} = $splicing_array[$j];
}


for (my $k=0; $k<scalar(@isoform_array); $k++ )
{
     my $isoform_array_string=$isoform_array[$k];
     $isoform_array_string=~s/^\s*|\s*$//g;    #like the String.trim() method
     my @string_isoform_array = split('\t', $isoform_array_string);
     my $trans_id = $string_isoform_array[0];
     my $gene_id = $string_isoform_array[1];
     my $gene_trans = $gene_id.":".$trans_id;
     $hash_gene_trans{$gene_trans} = $isoform_array[$k];
}

for(my $i=0; $i<scalar(@cufflinks_array); $i++ )  
{
  my $string_array=$cufflinks_array[$i];
  $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @string = split('\t', $string_array);
  my $attributes = $string[4]; 
  my $id_array = $string[3];
  $id_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @array_id= split(':', $id_array);
  my $gene_id = $array_id[0];
  my $tss_id = $array_id[1];
  my $trans_id = $array_id[2];
  my $gene_tss = $gene_id.":".$tss_id;  
  my $gene_trans = $gene_id.":".$trans_id;
  
  if((exists $hash_gene_tss{$gene_tss}) and (exists $hash_gene_trans{$gene_trans})){ 
        if($attributes eq $alternative_splicing_type){
            print $output_bed  $cufflinks_array[$i];
        }
    }
}
close $output_bed;

