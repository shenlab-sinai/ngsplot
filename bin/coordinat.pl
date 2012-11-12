#!/usr/bin/perl -w

##
#########################################
##
## File Name:
##   coordinat.pl
##
## Description:
##   parse the BED file and output the <attributes> characters's exonic coordinates
##
## Usage:
##   $0  <input.bed>  <attributes>  <output.bed>
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

# Program arguments.
if(@ARGV<3){
        print "Usage: $0 input.bed  attributes  output.bed\n";
        exit 0;
}

# Read gene structures into a hash table.
my $bed = $ARGV[0];
my $attributes_input = $ARGV[1];
open BED, "<$bed" or die "Open input bed file error:$!\n";

my $output = $ARGV[2];
my $output_bed;
open $output_bed, ">$output" or die "Open output bed file error: $!\n";


my %hash_tss = ();
my %hash_chrom = ();
my %hash_strand = ();

while(<BED>){
     chomp;
     my($chrom,$start,$end,$feature,$attributes,$strand) = split /\t/;
     if($attributes eq $attributes_input){
        my @feature_string = split(':', $feature);
        my $gene_id = $feature_string[0];
        my $tss_id = $feature_string[1];

        $gene_id=~s/^\s*|\s*$//g;    #like the String.trim() method
        $tss_id=~s/^\s*|\s*$//g;    #like the String.trim() method
        $chrom=~s/^\s*|\s*$//g;    #like the String.trim() method
        $strand=~s/^\s*|\s*$//g;    #like the String.trim() method

        my $gene_tss = $gene_id.":".$tss_id;
        $hash_chrom{$gene_tss} = $chrom;
        $hash_strand{$gene_tss} = $strand;

        if(exists $hash_tss{$gene_tss}){
            if($strand eq "+"){
               $hash_tss{$gene_tss}=$hash_tss{$gene_tss}.";".($start+1);
            }else{
               $hash_tss{$gene_tss}=$hash_tss{$gene_tss}.";".$end;        
            }
        }else{
            if($strand eq "+"){
               $hash_tss{$gene_tss}=$start+1;
            }else{
               $hash_tss{$gene_tss}=$end;
            }
        }
     }
}
close BED;

for my $gene_tss ( keys %hash_tss ) {
  my $chrom = $hash_chrom{$gene_tss};
  my $strand = $hash_strand{$gene_tss};
  my $tss_coordinat_string=$hash_tss{$gene_tss};
  my @tss_coordinat_array=split(';', $tss_coordinat_string);
  my $total_coordinat=0;
  for(my $i=0; $i<scalar(@tss_coordinat_array); $i++ ){
      $total_coordinat = $total_coordinat+$tss_coordinat_array[$i];
  }
  my $array_number = scalar(@tss_coordinat_array);
  my $real_tss_coordinat = int($total_coordinat/$array_number);
  my $bed_line=$chrom."\t".($real_tss_coordinat-1)."\t".$real_tss_coordinat."\t".$gene_tss."\t.\t".$strand."\n";
  print  $output_bed  $bed_line;
}

close $output_bed;
