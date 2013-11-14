#!/usr/bin/perl -w
use strict;

my $gemome = $ARGV[0];
my $region = $ARGV[1];

my $out_name = $ARGV[2];
my $out_avg_name = $ARGV[3];
my $out_hm_name = $ARGV[4];

my $database = $ARGV[5];

my $Randomly_sample=$ARGV[6];
my $GO=$ARGV[7];
my $CS= $ARGV[8];
my $FL= $ARGV[9];
my $MQ= $ARGV[10];
my $SE= $ARGV[11];
my $RB= $ARGV[12];
my $FC= $ARGV[13];
my $MW= $ARGV[14];
my $H=$ARGV[15];
my $flanking_region_option = $ARGV[16];

my $bam_file1 = $ARGV[17];
my $bam_file1_gene_list_detail = $ARGV[18];
my $image_title1 = $ARGV[19];
my $bam_file2 = $ARGV[20];
my $bam_file2_gene_list_detail = $ARGV[21];
my $image_title2 = $ARGV[22];

my $further_information = $ARGV[23];
my $interval_size = $ARGV[24];
my $flanking_region_size = $ARGV[25];

my $bam_file3 = $ARGV[26];
my $bam_file3_gene_list_detail = $ARGV[27];
my $image_title3 = $ARGV[28];

my $bam_file4 = $ARGV[29];
my $bam_file4_gene_list_detail = $ARGV[30];
my $image_title4 = $ARGV[31];

my $bam_file5 = $ARGV[32];
my $bam_file5_gene_list_detail = $ARGV[33];
my $image_title5 = $ARGV[34];

#my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;
#my $multiplot_conf_file="multiplot_conf".$out_base_name_png."txt";
open(FILE,">multiplot_conf_file.txt");

syswrite(FILE, "$bam_file1\t$bam_file1_gene_list_detail\t$image_title1\n");
syswrite(FILE, "$bam_file2\t$bam_file2_gene_list_detail\t$image_title2\n");
syswrite(FILE, "$bam_file3\t$bam_file3_gene_list_detail\t$image_title3\n");
syswrite(FILE, "$bam_file4\t$bam_file4_gene_list_detail\t$image_title4\n");
syswrite(FILE, "$bam_file5\t$bam_file5_gene_list_detail\t$image_title5\n");


close(FILE);

system("/home/galaxy/galaxy-dist/tools/multiplot/multiplot.genebody_exon_cgi.sh  $gemome $region   multiplot_conf_file.txt   $out_name  $out_avg_name  $out_hm_name  $database   $Randomly_sample  $GO  $CS $FL  $MQ  $SE  $RB  $FC  $MW  $H  $flanking_region_option  $further_information  $interval_size  $flanking_region_size");


