package Genome::Model::Tools::ViromeEvent::BlastN::PoolAndSplitSequence;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastN::PoolAndSplitSequence{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to pool together blast filtered files from human genomic blast stage and split up into files to NT blastN blast";
}

sub help_detail {
    "Tool to pool together blast filtered files from human genomic blast stage and split up into files to NT blastN blast";
}

sub execute {
    my $self = shift;
    unless( $self->pool_and_split_sequence( 'blastn', 500 ) ) {
	$self->error_message("Failed to successfully run pool and split sequence for blastn stage");
	return;
    }
    return 1;
}

1;
