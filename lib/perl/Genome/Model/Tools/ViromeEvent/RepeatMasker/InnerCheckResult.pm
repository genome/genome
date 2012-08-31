package Genome::Model::Tools::ViromeEvent::RepeatMasker::InnerCheckResult;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;
use Bio::SeqIO;

class Genome::Model::Tools::ViromeEvent::RepeatMasker::InnerCheckResult{
    is => 'Genome::Model::Tools::ViromeEvent',
    has =>[
        file_to_run => {
             is => 'String',  
            doc => 'files to rerun repeat masker', 
            is_input => 1,
        }
    ],
};

sub help_brief {
    return "module to make system call for repeat masker (decoupled for parallelization)";
}

sub help_detail {
    "Accepts a file for re-running repeat masker.  Assumes full path is included.";
}

sub execute {
    my $self = shift;

    my $file = $self->file_to_run;
    my $file_name = File::Basename::basename ($file);

    #ckeck if this is a re-attempt
    if (-s $file.'.masked') {
	my $input_count = 0;
	my $io = Bio::SeqIO->new(-format => 'fasta', -file => $file);
	while (my $seq = $io->next_seq) {
	    $input_count++;
	}
	my $masked_count = 0;
	my $masked_io = Bio::SeqIO->new(-format => 'fasta', -file => $file.'.masked');
	while (my $seq = $masked_io->next_seq) {
	    $masked_count++;
	}
	if($masked_count == $input_count) {
	    $self->log_event("Repeat masker already processed for $file_name");
	    return 1;
	}
    }

    $self->log_event("Running repeat masker on $file_name");

    #RUN REPEAT MASKER USING GZHAO'S LIBRARY
    my $com = "/gscmnt/sata835/info/medseq/virome/scripts_used_by_virome/RepeatMasker " . $file;

    if (system($com)) {
	$self->log_event("Repeat masker failed for $file_name");
	return;
    }
    $self->log_event("Finished running repeat masker on $file_name");
    return 1;
}

1;



