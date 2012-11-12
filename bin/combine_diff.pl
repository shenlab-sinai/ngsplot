#!/usr/bin/perl -w
##
#########################################
##
## File Name:
##   combine_diff.pl
##
## Description:
##   Combining the promoters.diff/splicing.diff file with dysregulated genes list file 
##   to get the both significant dysregulated genes list.
##
## Usage:
##   $0 <promoters.diff/splicing.diff> <dysregulated_list.file>  <co-significant_dysregulated_list.file> 
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
my $usage = "$0 <promoters.diff/splicing.diff> <dysregulated_list.file>  <co-significant_dysregulated_list.file> \n";

if(@ARGV<3){
        print "Usage: $0  <promoters.diff/splicing.diff> <dysregulated_list.file>  <co-significant_dysregulated_list.file>\n";
        exit 0;
}


my $infile = $ARGV[0];
my $dys_filename = $ARGV[1];
my $dys_co_significantname = $ARGV[2];

open(IN1, $infile) || die $usage ;
open(IN2, $dys_filename) || die $usage ;
my $dys_co_significant;

open $dys_co_significant, ">$dys_co_significantname" or die "Open output both significant dysregulated list file error: $!\n";

my @array=<IN1>;
my @dys_array=<IN2>;

for(my $i=0; $i<scalar(@array); $i++ )  
{
  my $string_array=$array[$i];
  $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @string = split('\t', $string_array);
  my $gene_id = $string[1];
    
  for (my $j=0; $j<scalar(@dys_array); $j++ )
  {
     my $dys_string=$dys_array[$j];
     $dys_string=~s/^\s*|\s*$//g;    #like the String.trim() method
     my @string_dys = split('\t', $dys_string);
     my $dys_gene_id = $string_dys[1];
     if($gene_id eq $dys_gene_id){ 
        print $dys_co_significant  $dys_array[$j];
     }
  }

}
close $dys_co_significant;

