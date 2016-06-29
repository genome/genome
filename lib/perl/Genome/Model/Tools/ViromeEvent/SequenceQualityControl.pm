package Genome::Model::Tools::ViromeEvent::SequenceQualityControl;

use strict;
use warnings;

use Genome;
use IO::File;
use Bio::SeqIO;
use File::Basename;

class Genome::Model::Tools::ViromeEvent::SequenceQualityControl{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    'filtering step prior to blast for virome pipeline';
}

sub help_detail {
    return <<"EOS"
This script will check each .masked file in the given directory.
Some sequences have only/lots of Ns because masked by RepeatMasker.
Remove sequences that do not have greater than 50 nt of consecutive
sequence without N.
EOS
}

sub execute {
    my $self = shift;

    my $dir = $self->dir;
    my $sample_name = basename ($dir);

    my $good_seq_file = $dir.'/'.$sample_name.'.fa.cdhit_out.masked.goodSeq';
    my $bad_seq_file  = $dir.'/'.$sample_name.'.fa.cdhit_out.masked.badSeq';

    my $good_io = Bio::SeqIO->new(-format => 'fasta', -file => ">$good_seq_file");
    my $bad_io  = Bio::SeqIO->new(-format => 'fasta', -file => ">$bad_seq_file");

    my $repeat_masker_dir = $dir.'/'.$sample_name.'.fa.cdhit_out_RepeatMasker';
    unless (-d $repeat_masker_dir) {
        $self->log_event("Failed to find repeat masker dir for sample: $sample_name");
        return;
    }

    my @reads_files = glob ("$repeat_masker_dir/$sample_name*.fa.masked");
    if ( not @reads_files ) {
        # no masked file means no reads were masked so get the original
        # input reads and proceed on
        @reads_files = glob("$repeat_masker_dir/$sample_name*.fa");
    }

    foreach my $file (@reads_files) {
	my $io = Bio::SeqIO->new(-format => 'fasta', -file => $file);
	while (my $seq = $io->next_seq) {
	    #LOOK FOR 50+ CONTIGUIOUS BASES(NONE-NS)
	    my $ge_50 = 0;
	    my @strings = split (/[n|x]/i, $seq->seq);
	    foreach my $string (@strings) {
		$ge_50 = 1 if (length $string >= 50);
	    }
	    if ($ge_50 == 1) {
		$good_io->write_seq($seq);
	    }
	    else {
		$bad_io->write_seq($seq);
	    }
	}
    }

    $self->log_event("Reads available for next stage for sample: $sample_name");
    return 1;
}

1;
