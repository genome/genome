package Genome::Model::Command::Define::DeNovoAssembly;

use strict;
use warnings;

use Genome;

require Carp;
use Regexp::Common;

class Genome::Model::Command::Define::DeNovoAssembly {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        center_name => {
            is => 'Text',
            valid_values => [qw/ WUGC LANL Baylor /],
            doc => 'Center name.'
        },
        import_location => {
            is => 'Text',
            is_optional => 1,
            doc => 'Directory to import assembly files from',
        },
    ]
};

sub type_specific_parameters_for_create {
    my $self = shift;
    my %p = (
        $self->SUPER::type_specific_parameters_for_create,
        center_name => $self->center_name,
    );
    $p{import_location} = $self->import_location if
        $self->import_location;
    #return ( center_name => $self->center_name );
    return %p;
}

sub listed_params {
    my $self = shift;
    return ($self->SUPER::listed_params, 'center_name', 'import_location');
}

1;

