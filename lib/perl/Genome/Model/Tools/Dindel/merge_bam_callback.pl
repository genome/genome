#!/usr/bin/env genome-perl

use strict;
use warnings;
use IO::File;
use Genome;
my ($bam_file, $output_prefix, $chromosome, $left_window_coord, $right_coord) = @ARGV;

my $output_file = $output_prefix . ".merged.sam";
my $merge_cmd = "samtools view $bam_file >> $output_file";

my $log_file = $output_prefix . ".windows.log";
my $log_fh = IO::File->new($log_file, ">>");
$log_fh->print("$chromosome\t$left_window_coord\t$right_coord\n");
$log_fh->close;

my $rv = Genome::Sys->shellcmd(cmd=>"$merge_cmd");

unlink($bam_file);


