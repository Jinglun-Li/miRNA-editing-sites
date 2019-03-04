#!/usr/bin/perl -w
use strict;
# This file removes the adepters and filters low quality reads
# Usage: $ perl process_reads.pl INPUT_Fastq_file OUTPUT_Fastq_file

### START OF INPUT
my $file_input=$ARGV[0];
my $file_output=$ARGV[1];
open FILE_R, "<", $file_input or die $!; #input fastq
open FILE_W, ">", $file_output or die $!; #output fastq (after filtering) 
my $minimal_quality=20; # the minimum quality score allowed in a given position 
my $max_number_of_strikes=3; # max number of bases with quality lower than $minimal_quality
my $error_min_final_length=15; # discard the read if it's too short after trimming of adaptors
my $error_max_final_length=28; # discard the read if it's too long after trimming of adaptors

my $must_be_at_start_of_adaptor3=5; # number of bases from the 3' adaptor that should appear in the read 
my $must_be_at_start_of_adaptor5=6; # changing this number require code changing
my $adaptor5="GTTCAGAGTTCTACAGTCCGACGATC"; # Illumina 5' adaptor 
my $adaptor3="ATCTCGTATGCCGTCTTCTGCTT"; # Illumina 3' adaptor; Can also be "TCGTATGCCGTCTTCTGCTT"
my $ascii_to_quality_convertor=33; # in recent Illumina CASAVA protocols (1.8+) this value in 33, in previous versions it was 64  

### END OF INPUT

my $number_of_seq=0;
my $seq_after_filter=0;
my $too_long_flag=0;
my $too_short_flag=0;
my $no_adaptor3_flag=0;
my $low_quality_flag=0;

my $string;
my @array;
my $value;
my $runner;
my $sequence;
my $ascii;
my $title_seq;
my $title_ascii;
my $include_flag;
my $start_at;
my $end_at;
my $i;
my $specific_length;
my $pattern;
my $temp_number;
my $strikes;
my $cmd;
my $length_identity_adaptor3;
my $temp_string;
my $currect_adaptor3;
my $possibly_miR;
my $possibly_barcode;
my $string_INFO;
my $temp_length;
my $start_time;
my $end_time;
my $temp_start_string;
my $temp_end_string;
my $end_adaptor5;
my $end_adaptor5_1;
my $end_adaptor5_2;
my $end_adaptor5_3;
my $end_adaptor5_4;
my $end_adaptor5_5;
my $miRNA_present;

$string=<FILE_R>;
$start_time=time();
while ($string)
{
    if ($string =~ /^@/)
    {
        $miRNA_present=0;
        $include_flag=1;
        $specific_length=0;
        $number_of_seq=$number_of_seq+1;
        $temp_number=$number_of_seq/10000;
        if ($temp_number !~ /\D/)
        {
            print "sequence number=",$number_of_seq,"\n";
        }
        $title_seq=$string;
        $string=<FILE_R>;
        $sequence=$string;
        
        # first screen
        $possibly_miR='';
        $end_adaptor5=substr($adaptor5,-$must_be_at_start_of_adaptor5);
        $end_adaptor5_1=substr($adaptor5,-1);
        $end_adaptor5_2=substr($adaptor5,-2);
        $end_adaptor5_3=substr($adaptor5,-3);
        $end_adaptor5_4=substr($adaptor5,-4);
        $end_adaptor5_5=substr($adaptor5,-5);        
        $temp_start_string=substr($adaptor3,0,$must_be_at_start_of_adaptor3);
#        if  ($sequence =~ /(\w+?)$temp_start_string/)
        if  ($sequence =~ /(\w+)/)
        {
            $possibly_miR=$1;
#            if ($possibly_miR =~ /$end_adaptor5(\w+)/)
#            {
#                $start_at=index($possibly_miR,$1);
#                $possibly_miR=$1;
#            }
#            elsif ($possibly_miR =~ /^$end_adaptor5_5(\w+)/)
#            {
#                $start_at=index($possibly_miR,$1);
#                $possibly_miR=$1;
#            }
#            elsif ($possibly_miR =~ /^$end_adaptor5_4(\w+)/)
#            {
#                $start_at=index($possibly_miR,$1);
#                $possibly_miR=$1;
#            }                
#            elsif ($possibly_miR =~ /^$end_adaptor5_3(\w+)/)
#            {
#                $start_at=index($possibly_miR,$1);
#                $possibly_miR=$1;
#            }
#            elsif ($possibly_miR =~ /^$end_adaptor5_2(\w+)/)
#            {
#                $start_at=index($possibly_miR,$1);
#                $possibly_miR=$1;
#            }
#            elsif ($possibly_miR =~ /^$end_adaptor5_1(\w+)/)
#            {
#                $start_at=index($possibly_miR,$1);
#                $possibly_miR=$1;
#            }                
#            else
#            {
                $start_at=0; 
#            }
            $temp_length=length($possibly_miR);
            if (($temp_length<=$error_max_final_length)&&($temp_length>=$error_min_final_length))
            {
                $end_at=$temp_length-1+$start_at; 
                $miRNA_present=1;
            }
            elsif ($temp_length>$error_max_final_length)
            {
                $too_long_flag=$too_long_flag+1;
                $include_flag=0;
            }
            elsif ($temp_length<$error_min_final_length)
            {
                $too_short_flag=$too_short_flag+1;
                $include_flag=0;
            }
       }
        
        if ($miRNA_present==1)
        {
            $specific_length=$end_at-$start_at+1;
            $string=<FILE_R>;            
        }
        else
        {
            if ($include_flag==1)
            {
                $include_flag=0;
            }
            $string=<FILE_R>;
            $string=<FILE_R>;
            $string=<FILE_R>;
        }
    }
    else
    {
        last;
    }
    if ($string =~ /^\+/)
    {
        $title_ascii=$string;
        $string=<FILE_R>;
        $ascii=substr($string,$start_at,$specific_length);
        @array=split(//,$string);
        $strikes=0;
        for($runner=$start_at;$runner<=$end_at;$runner++)
        {
            $value=ord($array[$runner])-$ascii_to_quality_convertor;
            if ($value<$minimal_quality)
            {
                $strikes=$strikes+1;           
            }
        }
        if ($strikes>$max_number_of_strikes)
        {
            $low_quality_flag=$low_quality_flag+1;
            $include_flag=0;
        }
        $string=<FILE_R>;
    }
    if ($include_flag)
    {
        $seq_after_filter=$seq_after_filter+1;
        print FILE_W $title_seq;
        print FILE_W $possibly_miR,"\n";
        print FILE_W $title_ascii;
        print FILE_W $ascii,"\n";
    }
}

close(FILE_R);
close(FILE_W);

$end_time=time();
$end_time=$end_time-$start_time;


print "number of sequences before filtering:\t",$number_of_seq,"\n";
print "number of sequences after filtering:\t",$seq_after_filter,"\n";
print "number of sequences filtered because they lack proper adaptors sequence:\t",$number_of_seq-$seq_after_filter-$too_long_flag-$too_short_flag-$low_quality_flag,"\n";
print "number of sequences filtered because they were too long:\t",$too_long_flag,"\n";
print "number of sequences filtered because they were too short:\t",$too_short_flag,"\n";
print "number of sequences filtered because they had low sequencing quality:\t",$low_quality_flag,"\n";
print "time taken to calc (sec):\t",$end_time,"\n";


