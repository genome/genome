package Genome::Model::Tools::ViromeEvent::RepeatMasker::SplitGivenNumberReads;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use Bio::SeqIO;
use File::Basename;

class Genome::Model::Tools::ViromeEvent::RepeatMasker::SplitGivenNumberReads{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    return "gzhao's Repeat Masker split given # of reads";
}

sub help_detail {
    return <<"EOS"
Given a fasta file, this script will split it to a number of files. Each 
file will contain given number of sequences. Generated files have the 
same name as the given file with numbered suffix .file0.fa .file1.fa ... 
etc All the generated files are placed in on subdirectory with the 
same name as the given file with "_RepeatMasker" suffix. 
EOS
}

sub execute
{
    my $self = shift;
    my $dir = $self->dir;
    my $sample_name = basename ($dir);

    #FIND CDHIT OUTPUT FILE
    my $cd_hit_result_file = $dir.'/'.$sample_name.'.fa.cdhit_out';
    unless (-s $cd_hit_result_file) {
	$self->log_event("Failed to find cd-hit result file for sample: $sample_name");
        #TODO - FAILUER HERE SHOULD CAUSE THE END OF SCREENING FOR THIS SAMPLE AND NOT RE-ATTEMPT
	return 0;
    }

    #SPLIT INTO 500 READS PER FILE
    my $output_dir = $dir.'/'.$sample_name.'.fa.cdhit_out_RepeatMasker';
    system ("mkdir $output_dir");
    unless (-d $output_dir) {
	$self->log_event("Failed to make repeat masker dir for sample: $sample_name");
	#TODO - FAILURE HERE SHOULD CAUSE THE END OF SCREENING FOR THIS SAMPLE AND NOT RE-ATTEMPT
	return 0;
    }

    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $cd_hit_result_file);
    unless ($in) {
	$self->log_event("Can not create Bio SeqIO for fasta file");
	#TODO - FAIURE HERE SHOULD CAUSE THE END OF SCREENING FOR THIS SAMPLE AND NOT RE-ATTEMPT
	return 0;
    }

    my $c = 0; my $n = 0; my $max = 2000;

    my $out_file = $output_dir.'/'.$sample_name.'.fa.cdhit_out_file'.$n.'.fa';

    my $out_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$out_file");
    unless (defined $out_io) {
	$self->log_event ("Can not create Bio SeqIO out for repeat masker splitting");
	return 0;
    }

    while (my $seq = $in->next_seq) {
	$c++;
	$out_io->write_seq($seq);
	if ($c == $max) { #REACHED MAX .. REDEFINE FILE NAME AND OUT IO HANDLE AND NOT RE-ATTEMPT
	    $c = 0;
	    $out_file = $output_dir.'/'.$sample_name.'.fa.cdhit_out_file'.++$n.'.fa';
	    $out_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$out_file");
	}
    }
    $self->log_event("Completed for sample: $sample_name");
    return 1;
}

1;
