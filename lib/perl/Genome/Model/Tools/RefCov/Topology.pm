package Genome::Model::Tools::RefCov::Topology;

use strict;
use warnings;

use Genome;

# ** NOTE ** Thu Jul  9 23:45:59 CDT 2009
# The coverage_topology() method will always return the revised version of the
# coverage-depth array as manipulated under the min_depth filter value.

# ** NOTE ** Fri Jul 10 00:34:31 CDT 2009
# Eventually, we will add plotting/graphing functions to this package for
# visualizing the coverage topology.

class Genome::Model::Tools::RefCov::Topology {
    has => [
        coverage => {
            is => 'ArrayRef',
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
    $self->coverage($coverage);
    return $self;
};

sub _set_depth_to_zero {
    my ($self, $position) = @_;
    $self->coverage->[$position] = 0;  # revise string
    return $self;
}


sub _revise_topology {
    my $self = shift;

    my $p = -1;
    POSITION:
    foreach my $position (@{ $self->coverage() }) {
        $p++;
        if ($position < $self->min_depth()) { $self->_set_depth_to_zero( $p ) }
    }

    return $self;
}


sub topology {
    my $self = shift;

    $self->_revise_topology();

    # Return revised coverage array reference.
    return $self->coverage();
}


sub save_topology {
    my ($self, $file) = @_;

    # Require a file path.
    if (!$file) { croak (__PACKAGE__ . ' save_topoloy requires a "file" argument.') }

    $self->_revise_topology();

    open (OUT, ">$file") or die 'could not open save file for topology';
    my $p = 0;
    foreach my $position (@{ $self->coverage() }) {
        $p++;
        print OUT join ("\t", $p, $position) . "\n";
    }
    close (OUT);

    return $self;
}


sub print_topology {
    my $self = shift;

    $self->_revise_topology();

    my $p = 0;
    foreach my $position (@{ $self->coverage() }) {
        $p++;
        print join ("\t", $p, $position) . "\n";
    }

    return $self;
}


1;  # End of package.
