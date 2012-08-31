package Genome::Model::Tools::RefCov::RelativeCoverage;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::RelativeCoverage {
    has => [
        coverage => {
            is => 'ArrayRef',
            is_optional => 1,
        },
        min_depth => {
            is => 'Integer',
            is_optional => 1,
            default_value => 0,
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $coverage = delete($params{coverage});
    my $self = $class->SUPER::create(%params);
    unless ($self) { return; }
    $self->coverage($coverage);
    return $self;
}

sub _set_depth_to_zero {
    my ($self, $position) = @_;
    $self->{_coverage}->[$position] = 0;  # revise string
    return $self;
}

sub relative_coverage {
    my $self = shift;
    $self->_revise_relative_coverage;
    return $self->{_relative_coverage};
}

sub _revise_relative_coverage {
    my $self = shift;

    my $positions = scalar(@{ $self->coverage() });
    my %relative_coverage;
    for (my $i = 0; $i < $positions; $i++) {
        my $depth = $self->coverage->[$i];
        unless ($depth) { next; }
        if ($self->min_depth && $depth < $self->min_depth) {
            $self->_set_depth_to_zero( $i )
        }
        my $relative_position = sprintf("%.02f",(($i + 1) / $positions));
        $relative_coverage{$relative_position} += $depth;
    }
    $self->{_relative_coverage} = \%relative_coverage;
}

sub print_relative_coverage {
    my $self = shift;

    my $relative_coverage = $self->relative_coverage;
    foreach my $relative_position (sort {$a <=> $b} keys %{$relative_coverage}) {
        print join ("\t", $relative_position, $relative_coverage->{$relative_position}) . "\n";
    }

    return $self;
}


1;  # End of package.
