#!/usr/bin/perl -w
##
#########################################
##
## File Name:
##   parse_diff.pl
##
## Description:
##   Parse the cuffdiff result (promoters.diff/splicing.diff)
##   and get the significant genes list.
##
## Usage:
##   $0 <promoters.diff/splicing.diff>  <significant_list.file> 
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
my $usage = "$0 <promoters.diff/splicing.diff>  <significant_list.file>\n";

if(@ARGV<2){
        print "Usage: $0 <promoters.diff/splicing.diff>  <significant_list.file>\n";
        exit 0;
}


my $infile = $ARGV[0];
my $significant_filename = $ARGV[1];


open(IN, $infile) || die $usage ;
my $significant_file;

open $significant_file, ">$significant_filename" or die "Open output significant list file error: $!\n";

my @array=<IN>;

for(my $i=1; $i<scalar(@array); $i++ )  #first line is the header
{
  my $string_array=$array[$i];
  $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @string = split('\t', $string_array);
  my $test_id = $string[0];
  my $gene_id = $string[1];
  my $gene_name = $string[2];
  my $significant = $string[13];
  if($significant eq "yes"){
        print $significant_file  $array[$i];
  }

}
close $significant_file;


