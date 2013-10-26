package Genome::Nomenclature;

use strict;
use warnings;

use Command::Dispatch::Shell;
use Genome;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use JSON::XS;

class Genome::Nomenclature {
    table_name => 'web.nomenclature',
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
            doc => 'Nomenclature name',
        },
        empty_equivalent => {
            is => 'Text',
            len => 25,
            is_optional => 1,
            column_name => 'default_value',
            doc => 'Empty-equivalent string (NA, n/a, etc)',
        },
        accepts_any_field => {
            is => 'Number',
            len => 1,
            default_value => 0,
            is_optional => 1,
            doc => 'Importer will create any field that doesnt exist in the nomenclature',
        },
        fields => {
            is => 'Genome::Nomenclature::Field',
            reverse_as => 'nomenclature',
            is_many => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'Nomenclatures',
};

sub create {
    my $class = shift;
    my %p = @_;

    if (exists $p{json}) {
        return $class->create_from_json($p{json});
    }

    $class->SUPER::create(@_);
}

sub create_from_json {
    my $class = shift;
    my $json = shift;

    my $nomenclature_raw = decode_json($json);

    my $nom = $class->create(
        name => $nomenclature_raw->{name},
        empty_equivalent => $nomenclature_raw->{'empty_equivalent'}
    );

    for my $rf (@{$nomenclature_raw->{fields}}) {
        my $f = Genome::Nomenclature::Field->create(name=>$rf->{name}, type=>$rf->{type}, nomenclature=>$nom);
        if ($rf->{type} eq 'enumerated') {
            for my $e (@{$rf->{enumerated_values}}) {
                my $enum = Genome::Nomenclature::Field::EnumValue->create(nomenclature_field=>$f, value=>$e);
            }
        }
    }
    return $nom;
}

sub json {
    my $self = shift;

    my $json = shift;
    if (!$json) {
        die "no JSON passed in";
    }
    my $nomenclature_raw = decode_json($json);
    if (!$nomenclature_raw) {
        die "no decodable JSON";
    }

    # step 1: update the name and empty_equivalent if necessary

    if ($nomenclature_raw->{name} ne $self->name) {
        $self->name($nomenclature_raw->{name});
    }

    if ( $nomenclature_raw->{'empty_equivalent'} ne $self->empty_equivalent() ) {
        $self->empty_equivalent($nomenclature_raw->{'empty_equivalent'});
    }

    # step 2: collect existing field ids so we know if we need to delete any fields.
    # if we see a field in the json then we'll delete it here, so at the end we'll have
    # a whole list of fields that need deleting altogether.
    my %fields_to_delete = map {$_->id, 1} Genome::Nomenclature::Field->get(nomenclature=>$self);


    for my $rf (@{$nomenclature_raw->{fields}}) {
        my $nom_field;
        if ($rf->{id}) {
            $nom_field = Genome::Nomenclature::Field->get($rf->{id});
        }

        # this must be an all new field.  let's add it and any enum values if needed.
        if (!$nom_field) {
            $nom_field = Genome::Nomenclature::Field->create(nomenclature=>$self, name=>$rf->{name}, type=>$rf->{type});
            for (@{$rf->{enumerated_values}}) {
                Genome::Nomenclature::Field::EnumValue->create(nomenclature_field=>$nom_field, value=>$_);
            }

            # nothing more to do since we've created the field and all that should go along with it.
            next;
        }


        # Is this an enumerated field?  Update the values for the field and remove unused ones.

        my %enum_records;
        warn $nom_field->type;
        if ($nom_field->type eq 'enumerated') {

            # as we scan the existing value ids we remove them, what's left in here needs to be deleted.
            my %enum_values_to_delete = map {$_->id,1} Genome::Nomenclature::Field::EnumValue->get(nomenclature_field=>$nom_field);

            my @value_ids = @{$rf->{enumerated_value_ids}};
            my @values = @{$rf->{enumerated_values}};
            for my $i (0...$#value_ids) {

                if ($value_ids[$i] == -1) {
                    Genome::Nomenclature::Field::EnumValue->create(nomenclature_field=>$nom_field, value=>$values[$i]);
                } else {
                    delete $enum_values_to_delete{$value_ids[$i]};
                    my $e = Genome::Nomenclature::Field::EnumValue->get($value_ids[$i]);
                    warn sprintf("Updating %s to %s", $e->value, $values[$i]);
                    $e->value($values[$i]);
                }
            }

            for(Genome::Nomenclature::Field::EnumValue->get(id=>[keys %enum_values_to_delete])) {
                $_->delete;
            }
        }

        if ($nom_field) {
            $nom_field->name($rf->{name});
            if ($rf->{type} eq 'enumerated' && $nom_field->type ne 'enumerated') {
                for my $e (@{$rf->{enumerated_values}}) {
                    my $enum = Genome::Nomenclature::Field::EnumValue->create(nomenclature_field=>$nom_field, value=>$e);
                }
            }
            if ($rf->{type} ne 'enumerated' && $nom_field->{type} eq 'enumerated') {
                for ($nom_field->enumerated_values) {
                    $_->delete;
                }
            }
            $nom_field->type($rf->{type})
        }

        delete $fields_to_delete{$nom_field->id};
    }

    for (Genome::Nomenclature::Field->get(id=>[keys %fields_to_delete])) {
        $_->delete;
    }

    return $self;
}

sub __display_name__ {
    shift->name;

}


1;
