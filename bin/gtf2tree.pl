#!/usr/bin/perl -w

##
#########################################
##
## File Name:
##   gtf2tree.pl
##
## Description:
##   Analyze a GTF file, annotate all exons and tss for each gene or transcript
##   and output exonic coordinates to a BED file.
##
## Usage:
##   $0  <input.gtf>   <output.bed>
##
## Author:
##   Dr. Xiaochuan Liu (xiaochuan.liu@mssm.edu)
##   Dr. Li Shen (li.shen@mssm.edu)
## Date:
##   Mon Dec 12 09:00:00 EST 2011
##
#########################################
#
#


use strict;

# Program arguments.
if(@ARGV<1){
        print "Usage: $0 input.gtf output.bed\n";
        exit 0;
}

# Read gene structures into a hash table.
my $gtf = $ARGV[0];
open GTF, "<$gtf" or die "Open input gtf file error:$!\n";
my %gene_table;
while(<GTF>){
    chomp;
    my($chrom,$source,$feature,$start,$end,$score,$strand,
                $frame,$attributes) = split /\t/;
    my %attr_table = ($attributes =~ /(\w+) \"(\S+)\"\;/g);
    my $gene_id = $attr_table{gene_id};
    my $tss_id = "";
    if(exists ($attr_table{tss_id})){$tss_id = $attr_table{tss_id};}
    my $transcript_id = $attr_table{transcript_id};
    if($tss_id ne "") {  
        if(exists $gene_table{$gene_id}){
             if(exists $gene_table{$gene_id}->{tss}{$tss_id}){
                   if(exists $gene_table{$gene_id}->{tss}{$tss_id}->{transcripts}{$transcript_id}){             
                         push @{$gene_table{$gene_id}->{tss}{$tss_id}->{transcripts}{$transcript_id}{exons}},
                                {'start' => $start, 'end' => $end, 'class' => 'canonical'};  
                    }else{
                        $gene_table{$gene_id}->{tss}{$tss_id}->{transcripts}{$transcript_id} = {
                                'exons' => [{'start' => $start, 'end' => $end, 'class' => 'canonical'}]
                        };            
                    }      
             }else{
                   $gene_table{$gene_id}->{tss}{$tss_id}={
                              'transcripts' => {
                                    $transcript_id => {
                                            'exons' => [{'start' => $start, 'end' => $end, 'class' => 'canonical'}]
                                }
                        }
                    }; 
              }
       }else{
                $gene_table{$gene_id} = {
                        'chrom' => $chrom,
                        'strand' => $strand,
                        'tss'  => {
                           $tss_id => {
                                'transcripts' => {
                                    $transcript_id => {
                                            'exons' => [{'start' => $start, 'end' => $end, 'class' => 'canonical'}]
                                }
                        }
                    }
                 }
             };  
        }
    }
}
close GTF;



# Annotate each exon and classify them as: promoter, polyA, variant, canonical and altBoundary.
while(my($gene_id, $gene_struct) = each %gene_table){
        my @tss_ids = keys %{$gene_struct->{tss}};
        my $strand = $gene_struct->{strand};
        if(@tss_ids > 1){  
            for my $i(0..$#tss_ids){
               my $tss_idA = $tss_ids[$i]; 
               my @trans_ids = keys %{$gene_struct->{tss}{$tss_idA}->{transcripts}};

               # Sort exons by start position
                foreach my $transcript(@trans_ids){
                        my @sorted_exons = sort {$a->{'start'} <=> $b->{'start'}} @{$gene_struct->{tss}{$tss_idA}->{transcripts}{$transcript}{exons}};
                        $gene_struct->{tss}{$tss_idA}->{transcripts}{$transcript}{exons} = [@sorted_exons];
                }

                
                   if(@trans_ids > 1){
                     for my $i(0..$#trans_ids-1){
                         for my $j(1..$#trans_ids){
                                my $idA = $trans_ids[$i]; 
                                my $idB = $trans_ids[$j];
                
                                &anno_exon($gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $gene_struct->{tss}{$tss_idA}->{transcripts}{$idB}{exons}, $strand);
                                &anno_exon($gene_struct->{tss}{$tss_idA}->{transcripts}{$idB}{exons}, $gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $strand);
                          }
                     }
                }else{
                my $idA = $trans_ids[0]; 
                &anno_exon($gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $strand); 
                }
            }   
        }else{
            my $tss_idA = $tss_ids[0];
            my @trans_ids = keys %{$gene_struct->{tss}{$tss_idA}->{transcripts}};
            # Sort exons by start position
            foreach my $transcript(@trans_ids){
                my @sorted_exons = sort {$a->{'start'} <=> $b->{'start'}} @{$gene_struct->{tss}{$tss_idA}->{transcripts}{$transcript}{exons}};
                $gene_struct->{tss}{$tss_idA}->{transcripts}{$transcript}{exons} = [@sorted_exons];
            }

                if(@trans_ids > 1){
                     for my $i(0..$#trans_ids-1){
                         for my $j(($i+1)..$#trans_ids){
                                my $idA = $trans_ids[$i];
                                my $idB = $trans_ids[$j];
                                &anno_exon($gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $gene_struct->{tss}{$tss_idA}->{transcripts}{$idB}{exons}, $strand);
                                &anno_exon($gene_struct->{tss}{$tss_idA}->{transcripts}{$idB}{exons}, $gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $strand);
                          }
                     }
                }else{
                my $idA = $trans_ids[0];
                &anno_exon($gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $gene_struct->{tss}{$tss_idA}->{transcripts}{$idA}{exons}, $strand);
                }
        }            
}



# Output exonic coordinates with annotation.
my $bed = $ARGV[1];
my $beh;
open $beh, ">$bed" or die "Open output BED file error: $!\n";
while(my($gene_id, $gene_struct) = each %gene_table){
         &output_gene($gene_id, $gene_struct, $beh);
}
close $beh;


###### Start Subroutines #######
# Function for exon annotation.
sub anno_exon{
        my($exons_ref,$exons_working,$strand) = @_;
        my $ref_exon_n = @{$exons_ref};
        my $working_exon_n = @{$exons_working}; 
        my $w = 0;      # counter for working exon.
        ## Function to determine whether two exons overlap.
        sub overlapExon{
                my($e1,$e2) = @_;
                return ($e1->{start} <= $e2->{end} and $e1->{end} >= $e2->{start});
        }
        foreach my $working_exon(@{$exons_working}){
                $w++; 
                # The first and last working exons are directly assigned class.
                if(($w==1 and $strand eq '+') or ($w==$working_exon_n and $strand eq '-')){
                        $working_exon->{class} = 'promoter';   
                }elsif(($w==1 and $strand eq '-') or ($w==$working_exon_n and $strand eq '+')){
                        $working_exon->{class} = 'polyA';
                }else{  # deal with working exon in the middle.
                        next if $working_exon->{class} eq 'variant';    # "variant" has priority over other types.
                        my $r = 0;      # iterator for reference exon.
                        # Mutually exclusive conditions...
                        while($r < $ref_exon_n and $working_exon->{end} >= $exons_ref->[$r]{start}){
                                $r++;
                        }
                        next if $r==0;  # working exon is before/after the reference promoter/polyA.
                        next if ($r==$ref_exon_n and $working_exon->{start} > $exons_ref->[$r-1]{end}); # similar as above.
                        if(!overlapExon($working_exon, $exons_ref->[$r-1])){    # variant exon.
                                $working_exon->{class} = 'variant';
                        }else{  # overlap: alternative boundaries.
                                next if $working_exon->{class} eq 'altBoth';    # "altBoth" has 2nd priority after "variant".
                                my $tentative_class;    # classification based on current overlapping info.
                                if($r==1 and $working_exon->{end} != $exons_ref->[$r-1]{end}){  # overlapping left-most ref exon.
                                        $tentative_class = $strand eq '+'? 'altDonor' : 'altAcceptor';
                                }elsif($r==$ref_exon_n and $working_exon->{start} != $exons_ref->[$r-1]{start}){        # overlapping right-most ref exon.
                                        $tentative_class = $strand eq '+'? 'altAcceptor' : 'altDonor';
                                }elsif($r > 1 and $r < $ref_exon_n){    # overlapping middle ref exon.
                                        if($working_exon->{start} != $exons_ref->[$r-1]{start} and
                                                $working_exon->{end} == $exons_ref->[$r-1]{end}){
                                                $tentative_class = $strand eq '+'? 'altAcceptor' : 'altDonor';
                                        }elsif($working_exon->{start} == $exons_ref->[$r-1]{start} and
                                                $working_exon->{end} != $exons_ref->[$r-1]{end}){
                                                $tentative_class = $strand eq '+'? 'altDonor' : 'altAcceptor';
                                        }elsif($working_exon->{start} != $exons_ref->[$r-1]{start} and
                                                $working_exon->{end} != $exons_ref->[$r-1]{end}){
                                                $tentative_class = 'altBoth';
                                        }
                                }
                                # Must consider previous class and current assignment.
                                if(defined $tentative_class){
                                        if(($working_exon->{class} eq 'altDonor' and $tentative_class eq 'altAcceptor') or
                                                ($working_exon->{class} eq 'altAcceptor' and $tentative_class eq 'altDonor')){
                                                $working_exon->{class} = 'altBoth';     # "upgrade" to both alternatives.
                                        }else{
                                                $working_exon->{class} = $tentative_class;      # tentative class prevails.
                                        }
                                }
                        }
                }
        }
}

# Output a whole gene structure.
sub output_gene{
        my($gene_id,$gene_struct,$beh) = @_;
        while(my($tss_id,$tss_struct) = each %{$gene_struct->{tss}}){
                while(my($transcript_id,$trans_struct) = each %{$gene_struct->{tss}{$tss_id}->{transcripts}}){
                     my $trans_desc = join(':', $gene_id, $tss_id, $transcript_id);
                        &output_transcript($gene_struct->{chrom}, $trans_desc,$gene_struct->{strand}, $trans_struct->{exons},$beh);
                }
        } 
}

# Output a transcript structure.
# Note: Convert 1-based coordinates in GTF to 0-based in BED!
sub output_transcript{
        my($chrom,$trans_desc,$strand,$exons,$h) = @_;
        foreach my $e(@{$exons}){
                # Output exon.
                print $h join("\t", $chrom,$e->{start}-1,$e->{end},$trans_desc,
                        'exon:' . $e->{class},$strand), "\n";

        }
}





