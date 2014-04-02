package Genome::Test::Factory::InstrumentData::MergedAlignmentResult;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;

#our @required_params = qw(library_id);

sub generate_obj {
    my $self = shift;
    return Genome::InstrumentData::AlignmentResult::Merged->__define__(@_);
}


1;
