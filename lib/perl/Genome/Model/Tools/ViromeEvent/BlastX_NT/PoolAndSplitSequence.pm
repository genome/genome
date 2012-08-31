package Genome::Model::Tools::ViromeEvent::BlastX_NT::PoolAndSplitSequence;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastX_NT::PoolAndSplitSequence{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to pool together blast filtered files from NT blastN stage and split up into files for NT blastX blast";
}

sub help_detail {
    "Tool to pool together blast filtered files from NT blastN stage and split up into files for NT blastX blast";
}

sub execute {
    my $self = shift;
    unless( $self->pool_and_split_sequence( 'blastx_nt', 250 ) ) {
	$self->error_message("Failed to successfully run pool and split sequence for blastx_nt stage");
	return;
    }
    return 1;
}

1;
