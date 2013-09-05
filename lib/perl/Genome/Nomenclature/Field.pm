package Genome::Nomenclature::Field;

use strict;
use warnings;

use Command::Dispatch::Shell;
use Genome;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use JSON::XS;

class Genome::Nomenclature::Field {
    table_name => 'web.nomenclature_field',
    id_by => [
        id => {
            is => 'Text',
            len => 255,
        },
    ],
    has => [
        name => {
            is => 'Text',
            len => 255,
            doc => 'Nomenclature field name',
        },
        type => {
            is => 'Text',
            len => 255,
            doc => 'Nomenclature field type',
        },
        nomenclature_id => {
            is => 'Text',
            len => 255,
        },
        nomenclature => {
            is => 'Genome::Nomenclature',
            id_by => 'nomenclature_id',
            constraint_name => 'GNF_FK',
        },
        nomenclature_name => {
            via => 'nomenclature',
            to => 'name',
        },
        enumerated_values => {
            is => 'Genome::Nomenclature::Field::EnumValue',
            reverse_as => 'nomenclature_field',
            is_many => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'Nomenclature::Fields',
};

sub __display_name__ {
    my $self = shift;
    sprintf("%s (%s)", $self->name, $self->type);
}

sub delete {
    my $self = shift;

    for ($self->enumerated_values) {
        $_->delete;
    }
    $self->SUPER::delete(@_);
}


1;
