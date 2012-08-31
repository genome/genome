package Genome::Model::Tools::ViromeEvent::BlastX_Viral::PoolAndSplitSequence;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastX_Viral::PoolAndSplitSequence{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to pool together NT blastX filtered files and split out into files for viral blastX";
}

sub help_detail {
    "Tool to pool together NT blastX filtered files and split out into files for viral blastX";
}

sub execute {
    my $self = shift;
    unless( $self->pool_and_split_sequence( 'blastx_viral', 500 ) ) {
	$self->error_message("Failed to successfully run pool and split sequence for blastx_viral stage");
	return;
    }
    return 1;
}

1;
