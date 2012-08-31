package Genome::Model::Tools::ViromeEvent::SequenceQualityControl;

use strict;
use warnings;

use Genome;
use Workflow;
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

    my @masked_files = glob ("$repeat_masker_dir/$sample_name*.fa.masked");
    unless (scalar @masked_files > 0) {
	#CHECK FOR *cat.all FILE .. IF PRESENT IT MEANS ALL READS WERE PROCESSED OUT
	my @cat_files = glob("$repeat_masker_dir/$sample_name*.cat.all");
	if (@cat_files > 0) {
	    $self->log_event("All reads have been filtered out for $sample_name");
	    return 1;
	}
	#OTHERWISE IT FAILED
        $self->log_event("No masked fastas for repeat masker run found for sample: $sample_name");
        return;
    }

    foreach my $file (@masked_files) {
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

    $self->log_event("Repeat masker Seq qual control completed for sample: $sample_name");
    return 1;
}

1;
