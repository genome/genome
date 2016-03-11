package Genome::InstrumentData::Command::Import::WorkFlow::Builder::Bam;

use strict;
use warnings;

use Genome;

require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::Builder::Bam {
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::Builder',
};

sub _steps_to_build_workflow {
    my $self = shift;

    my @steps = ( 'verify not imported', 'sort bam', );
    push @steps, 'downsample bam' if $self->work_flow_inputs->instrument_data_properties->{downsample_ratio};
    push @steps, 'sanitize and split bam', 'create instrument data';

    return @steps;
}

1;

