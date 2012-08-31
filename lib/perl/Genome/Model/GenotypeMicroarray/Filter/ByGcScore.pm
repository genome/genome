package Genome::Model::GenotypeMicroarray::Filter::ByGcScore;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Filter::ByGcScore {
    has => [
        min => {
            is => 'Text',
        },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    if ( not defined $self->min ) {
        $self->error_message("No min given to filter by gc score.");
        $self->delete;
        return;
    }

    if ( $self->min < 0 ) {
        $self->error_message("Invalid min (".$self->min.") given to filter by gc score.");
        $self->delete;
        return;
    }

    return $self;
}

sub filter { # does this variant PASS the filter?
    my ($self, $variant) = @_;

    return 1 if not defined $variant->{gc_score}; # do not filter if not defiend
    return 1 if $variant->{gc_score} == -1; # do not filter if -1 => this means no score
    return 1 if $variant->{gc_score} >= $self->min; # do not filter if greater than the min

    return; # filter!
}

1;

