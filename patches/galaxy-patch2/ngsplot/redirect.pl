#!/usr/bin/perl -w
use strict;
my $usage = "$0 <in1>\n";

my $infile1 = $ARGV[0];
open(IN1, $infile1) || die $usage ;

my @array1=<IN1>;
#print STDERR $infile1;
for(my $i=0; $i<scalar(@array1); $i++ )
{
  my $string=$array1[$i];
  print STDERR $string;
}
my @inf = split('\/', $infile1);
my $length=scalar(@inf);
print STDERR $inf[$length-1]."\n";

