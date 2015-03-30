#!/usr/bin/perl -w
use strict;
use File::Basename 'dirname';
use File::Spec;
use Cwd 'abs_path';

my @inputs = @ARGV;
my @inputs2 = @inputs;

my $genome_name = shift(@inputs);
my $genomic_region_source_type__genomic_region = shift(@inputs);
my $genomic_region_source_type__further_information = shift(@inputs);
my $genomic_region_source_type__interval_size = shift(@inputs);
my $genomic_region_source_type__flanking_region_option_source_type__flanking_region_option = shift(@inputs);
my $genomic_region_source_type__flanking_region_option_source_type__flanking_region_size = shift(@inputs);
my $numsamples = shift(@inputs);
my $usepairs = shift(@inputs);

#print STDERR "inputs @inputs\n";

my $randfile = rand(100)."\.config\.txt";
my $randfile2 = $randfile;
$randfile2 =~ s/config\.txt/logfile/gm;

my $outfile = File::Spec->catfile(abs_path(dirname(__FILE__)),"$randfile");
open(FILE,">$outfile");

for (my $i=1;$i<=$numsamples;$i++) {
   my $bamfile=shift(@inputs);
   my $reffile=shift(@inputs);
   my $usegenelist=shift(@inputs);
   my $genelist=shift(@inputs);
   my $title=shift(@inputs);
   my $fraglen=shift(@inputs);
   my $color=shift(@inputs);   
   if ($usepairs eq 'yes') {
      syswrite(FILE, "$bamfile\:$reffile\t$genelist\t$title\t$fraglen\t$color\n");
   }else {
      syswrite(FILE, "$bamfile\t$genelist\t$title\t$fraglen\t$color\n");
   }
}
close(FILE);

my $gene_database = shift(@inputs);
my $randomly_sample = shift(@inputs);
my $gene_order = shift(@inputs);
my $knc = shift(@inputs);
my $mit = shift(@inputs);
my $nrs = shift(@inputs);
my $chunk_size = shift(@inputs);
my $quality_requirement = shift(@inputs);
my $standard_error = shift(@inputs);
my $radius_size = shift(@inputs);
my $flooding_fraction = shift(@inputs);
my $smooth_method = shift(@inputs);
my $shaded_area = shift(@inputs);
my $out_name = shift(@inputs);
my $out_avg_name = shift(@inputs);
my $out_hm_name = shift(@inputs);


my $G = $genome_name;
my $R = $genomic_region_source_type__genomic_region;
my $C = $outfile;
my $O = $out_name;
my $O2 = $out_avg_name;
my $O3 = $out_hm_name;
my $F = $genomic_region_source_type__further_information;
my $D = $gene_database;
my $L = $genomic_region_source_type__flanking_region_option_source_type__flanking_region_size;
my $N = $genomic_region_source_type__flanking_region_option_source_type__flanking_region_size;
my $RB = $radius_size;
my $S = $randomly_sample;
my $CS = $chunk_size;
my $MQ = $quality_requirement;
my $IN = $genomic_region_source_type__interval_size;
my $SE = $standard_error;
my $MW = $smooth_method;
my $H = $shaded_area;
my $GO = $gene_order;
my $KNC = $knc;
my $MIT = $mit;
my $NRS = $nrs;
my $FC = $flooding_fraction;

if ($GO eq 'km') {
   $GO = "$GO -KNC $KNC -MIT $MIT -NRS $NRS";
}

my $logfile = File::Spec->catfile(abs_path(dirname(__FILE__)),"$randfile2");
open(FILE2,">>$logfile");
my $cmd5="pwd >> $logfile 2>&1";
system($cmd5);

my $cmd='';
if (($R eq 'tss')||($R eq 'tes')) {
   $cmd = "ngs.plot.r -Galaxy 1 -P 0 -G $G -R $R -C $C -O $O -O2 $O2 -O3 $O3 -D $D -L $L -S $S -GO $GO -CS $CS -MQ $MQ -SE $SE -RB $RB -FC $FC -MW $MW -H $H >> $logfile 2>&1";
}elsif ($genomic_region_source_type__flanking_region_option_source_type__flanking_region_option eq 'flanking_region_size') {
   if ($IN eq 'automatic') {
      $cmd = "ngs.plot.r -Galaxy 1 -P 0 -G $G -R $R -C $C -O $O -O2 $O2 -O3 $O3 -D $D -L $L -S $S -GO $GO -CS $CS -MQ $MQ -SE $SE -RB $RB -FC $FC -MW $MW -H $H >> $logfile 2>&1";
   }else {
      $cmd = "ngs.plot.r -Galaxy 1 -P 0 -G $G -R $R -C $C -O $O -O2 $O2 -O3 $O3 -D $D -L $L -S $S -GO $GO -CS $CS -MQ $MQ -SE $SE -RB $RB -FC $FC -MW $MW -H $H -IN $IN  >> $logfile 2>&1";
   }
}elsif ($genomic_region_source_type__flanking_region_option_source_type__flanking_region_option eq 'flanking_floating_size') {
   if ($IN eq 'automatic') {
      $cmd = "ngs.plot.r -Galaxy 1 -P 0 -G $G -R $R -C $C -O $O -O2 $O2 -O3 $O3 -D $D -N $N -S $S -GO $GO -CS $CS -MQ $MQ -SE $SE -RB $RB -FC $FC -MW $MW -H $H  >> $logfile 2>&1";
   }else {
      $cmd = "ngs.plot.r -Galaxy 1 -P 0 -G $G -R $R -C $C -O $O -O2 $O2 -O3 $O3 -D $D -N $N -S $S -GO $GO -CS $CS -MQ $MQ -SE $SE -RB $RB -FC $FC -MW $MW -H $H -IN $IN  >> $logfile 2>&1";
   }
}
my $cmd2="cp data.zip $O";
my $cmd3="rm $outfile";
my $cmd4="rm $logfile";
syswrite(FILE2, "\n$cmd\n");
syswrite(FILE2, "\n$cmd2\n");
syswrite(FILE2, "\n$cmd3\n");
syswrite(FILE2, "\n$cmd4\n\n");

system($cmd);
system($cmd2);
system($cmd3);
system($cmd4);

close(FILE2);


