#!/usr/bin/perl -w
use strict;
use File::Basename 'dirname';
use File::Spec;
use Cwd 'abs_path';

my $input_file = $ARGV[0];
my $output_type = $ARGV[1];
my $output_filename = $ARGV[2];

my $cmd1 = "cp $input_file data.zip";
my $cmd2 = "";
my $cmd3 = "mv output.pdf $output_filename";

my $FS = "";
if ($output_type eq 'prof') {
   my $WD = "-WD $ARGV[3]";
   my $HG = "-HG $ARGV[4]";
   my $SE = "-SE $ARGV[5]";
   my $H = "-H $ARGV[6]";
   my $MW = "-MW $ARGV[7]";
   my $Y = "-Y $ARGV[8]";
   my $YAS = "-YAS $ARGV[9]";
   my $LEG = "-LEG $ARGV[10]";
   my $BOX = "-BOX $ARGV[11]";
   my $VLN = "-VLN $ARGV[12]";
   my $XYL = "-XYL $ARGV[13]";
   my $LWD = "-LWD $ARGV[14]";
   $FS = "-FS $ARGV[15]";

   if ($Y =~ /no$/gm) {
      $cmd2 = "replot.r $output_type -I data.zip -O output $WD $HG $SE $H $MW $LEG $BOX $VLN $XYL $LWD $FS";
   }else {
      $cmd2 = "replot.r $output_type -I data.zip -O output $WD $HG $SE $H $MW $LEG $BOX $VLN $XYL $LWD $FS $YAS";
   }
}elsif ($output_type eq 'heatmap') {
   my $GO = "-GO $ARGV[3]";
   my $KNC = "-KNC $ARGV[4]";
   my $MIT = "-MIT $ARGV[5]";
   my $NRS = "-NRS $ARGV[6]";
   my $LOW = "-LOW $ARGV[7]";
   my $RR = "-RR $ARGV[8]";
   my $FC = "-FC $ARGV[9]";
   my $CO = "-CO $ARGV[10]";
   my $CD = "-CD $ARGV[11]";
   $FS = "-FS $ARGV[12]";
   my $DUM1 = $ARGV[13];
   my $DUM2 = $ARGV[14];
   my $DUM3 = $ARGV[15];

   if ($CO =~ /default/gm) { $CO = ""; }

   if ($GO =~ /km$/gm) {
      $cmd2 = "replot.r $output_type -I data.zip -O output $LOW $RR $FC $CO $CD $FS $GO $KNC $MIT $NRS"; 
   }else {
      $cmd2 = "replot.r $output_type -I data.zip -O output $LOW $RR $FC $CO $CD $FS $GO";
   }
}
system($cmd1);
system($cmd2);
system($cmd3);


