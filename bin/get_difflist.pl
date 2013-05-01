#!/usr/bin/perl -w

##
#########################################
##
## File Name:
##   get_diff.pl
##
## Description:
##   Parse the cuffdiff result and get the up-regulated and down-regulated genes list.
##
## Usage:
##   $0 <cuffdiff_result.file> <Treat:Control/Control:Treat>  <FoldChange>  <up-regulated_list.file>  <down-regulated_list.file>
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
my $usage = "$0 <cuffdiff_result.file> <condition>  <FoldChange>  <up-regulated_list.file>   <down-regulated_list.file>\n\nCondition:\n\t1) Treat:Control\n\t2) Control:Treat\n";

if(@ARGV<4){
        print "Usage: $0 <cuffdiff_result.file> <condition>  <FoldChange>  <up-regulated_list.file>  <down-regulated_list.file>\n\nCondition:\n\t1) Treat:Control\n\t2) Control:Treat\n";
        exit 0;
}


my $infile = $ARGV[0];
my $condition = $ARGV[1];
my $FC = $ARGV[2];
my $up_filename = $ARGV[3];
my $down_filename = $ARGV[4];

if(!(($condition eq "Treat:Control")or($condition eq "Control:Treat")))
{
   print "Usage: $0 <cuffdiff_result.file> <condition>  <FoldChange> <up-regulated_list.file>  <down-regulated_list.file>\n\nCondition:\n\t1) Treat:Control\n\t2) Control:Treat\n";
        exit 0;
}


if(!($FC > 0))
{
    print "<FoldChange> must be a positive integer !\n";
    print "Usage: $0 <cuffdiff_result.file> <condition>  <FoldChange> <up-regulated_list.file>  <down-regulated_list.file>\n\nCondition:\n\t1) Treat:Control\n\t2) Control:Treat\n";
        exit 0;
}

my $base=2;
my $log_FC = log($FC)/log($base); 


open(IN, $infile) || die $usage ;
my $up_file;
my $down_file;
open $up_file, ">$up_filename" or die "Open output up-regulated_list file error: $!\n";
open $down_file, ">$down_filename" or die "Open output down-regulated_list file error: $!\n";

my @array=<IN>;

for(my $i=1; $i<scalar(@array); $i++ )  #first line is the header
{
  my $string_array=$array[$i];
  $string_array=~s/^\s*|\s*$//g;    #like the String.trim() method
  my @string = split('\t', $string_array);
  my $test_id = $string[0];
  my $gene_id = $string[1];
  my $gene_name = $string[2];
  my $log_value = $string[9];
  my $p_value = $string[12];
  my $significant = $string[13];
  if($p_value < 0.05){
#  if($significant eq "yes"){
     if($condition eq "Control:Treat"){
        if($log_value > $log_FC ){ 
           print $up_file  $array[$i];
        }  
        if($log_value < (($log_FC)*(-1))){
           print $down_file  $array[$i];
        }
     }
     if($condition eq "Treat:Control"){
        if($log_value <  (($log_FC)*(-1))){
           print $up_file  $array[$i];
        }
        if($log_value >  $log_FC){
           print $down_file  $array[$i];
        }
     }
  }

}
close $up_file;
close $down_file;

