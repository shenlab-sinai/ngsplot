#!/usr/bin/perl -w
##
#########################################
##
## File Name:
##   difflist2bed.pl
##
## Description:
##   Parse the bed file and get the bed annotation for dysregulated genes list.
##
## Usage:
##   $0 <cufflinks_bed_file> <dysregulated_list.file>  <dysregulated_list.bed> 
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
my $usage = "$0 <cufflinks_bed_file> <dysregulated_list.file>  <dysregulated_list.bed>  \n";

if(@ARGV<3){
        print "Usage: $0 <cufflinks_bed_file> <dysregulated_list.file> <dysregulated_list.bed>\n";
        exit 0;
}


my $infile = $ARGV[0];
my $dys_filename = $ARGV[1];
my $dys_bedname = $ARGV[2];

open(IN1, $infile) || die $usage ;
open(IN2, $dys_filename) || die $usage ;
my $dys_bed;
my %hash_significant = ();

open $dys_bed, ">$dys_bedname" or die "Open output dysregulated list bed file error: $!\n";

my @array=<IN1>;
my @dys_array=<IN2>;


for (my $j=0; $j<scalar(@dys_array); $j++ )
{
     my $dys_string=$dys_array[$j];
     $dys_string=~s/^\s*|\s*$//g;    #like the String.trim() method
     my @string_dys = split('\t', $dys_string);
     my $dys_tss_id = $string_dys[0];
     my $dys_gene_id = $string_dys[1];
     my $gene_tss = $dys_gene_id.":".$dys_tss_id;
     
     $hash_significant{$gene_tss} = $dys_array[$j];
}


for(my $i=0; $i<scalar(@array); $i++ )  
{
  my $string_array=$array[$i];
  $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @string = split('\t', $string_array);
  my $id_array = $string[3];
  $id_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @array_id= split(':', $id_array);
  my $gene_id = $array_id[0];
  my $tss_id = $array_id[1];
  my $gene_tss = $gene_id.":".$tss_id;  
  
  if(exists $hash_significant{$gene_tss}){ 
        print $dys_bed  $array[$i];
    }
}
close $dys_bed;

