package Genome::Model::Tools::ViromeEvent::BlastHumanGenome::SplitGivenNumberReads;

use strict;
use warnings;

use Genome;
use IO::File;

use Bio::SeqIO;

class Genome::Model::Tools::ViromeEvent::BlastHumanGenome::SplitGivenNumberReads{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to pool blast filtered files from previous stage and create a new set of files for blast";
}

sub help_detail {
    "Tool to pool blast filtered files from previous stage and create a new set of files for blast";
}

sub execute {
    my $self = shift;
    my $dir = $self->dir;
    my $sample_name = File::Basename::basename ($dir);

   #DIRECTORY TO PUT SPLIT FILES INTO
    my $output_dir = $dir.'/'.$sample_name.'.fa.cdhit_out.masked.goodSeq_HGblast';
    Genome::Sys->create_directory( $output_dir ) unless -d $output_dir;

    #FILE TO SPLIT
    my $good_seq_file = $dir.'/'.$sample_name.'.fa.cdhit_out.masked.goodSeq';
    #IF FILE DOESN'T EXIST SOMETHING WENT WRONG
    unless (-e $good_seq_file) {
	$self->log_event("Failed to find repeat masker good seq file for $sample_name");
	return;
    }
    #IF FILE IS ZERO SIZE ALL READS HAVE BEEN PROCESSED OUT
    if (-s $good_seq_file == 0) {
	$self->log_event("No reads available for further processing for $sample_name");
	return 1;
    }

    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $good_seq_file);
    unless ($in) {
	$self->log_event("Failed to create Bio::SeqIO for human blast");
	return;
    }

    my $c = 0; my $n = 0; my $max = 500;

    my $out_file = $output_dir.'/'.$sample_name.'.fa.cdhit_out.masked.goodSeq_file'.$n.'.fa';

    my $out_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$out_file");
    unless (defined $out_io) {
	$self->log_event("Failed to create Bio SeqIO for HG blast splitting");
	return;
    }

    while (my $seq = $in->next_seq) {
	$c++;
	$out_io->write_seq($seq);
	if ($c == $max) {
	    $c = 0;
	    $out_file = $output_dir.'/'.$sample_name.'.fa.cdhit_out.masked.goodSeq_file'.++$n.'.fa';
	    $out_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$out_file");
	}
    }

    unless (-s $out_file) {
	$self->log_event("Failed to create HG blastN input file or no valid data available");
	return;
    }

    $self->log_event("Split reads for HG blast completed for $sample_name");

    return 1;
}

1;
