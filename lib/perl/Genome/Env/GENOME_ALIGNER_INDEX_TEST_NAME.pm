package Genome::Env::GENOME_ALIGNER_INDEX_TEST_NAME;

use strict;
use warnings;

use Genome;

our $VERSION = $Genome::VERSION;

=pod

=head1 NAME

GENOME_ALIGNER_INDEX_TEST_NAME

=head1 DESCRIPTION

The GENOME_ALIGNER_INDEX_TEST_NAME environment variable provides for making
a new AlignerIndex with what would otherwise be identical parameters to an
existing result, in order to run tests or to "move" results out of the way to be
regenerated.

=head1 DEFAULT VALUE

 undef

=cut

1;
