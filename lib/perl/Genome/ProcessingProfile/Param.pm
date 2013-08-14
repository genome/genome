package Genome::ProcessingProfile::Param;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::Param {
    table_name => 'model.processing_profile_param',
    type_name => 'processing profile param',
    id_by => [
        processing_profile_id => {
            is => 'Text',
            len => 32,
            constraint_name => 'PPP_PP_FK',
        },
        name => {
            is => 'Text',
            len => 255,
            column_name => 'PARAM_NAME',
        },
        value_id => {
            is => 'Text',
            len => 1000,
            column_name => 'PARAM_VALUE',
        },
    ],
    has => [
        value_class_name => {
            is => 'Text',
            len => 255,
        },
        value_obj => {
            is => 'UR::Object',
            id_by => 'value_id',
            id_class_by => 'value_class_name',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
        },
    ],
    has_deprecated_optional => [
        deprecated_name => {
            is => 'Text',
            len => 255,
            column_name => 'NAME',
        },
        deprecated_value_class_id => {
            is => 'Text',
            len => 1000,
            column_name => 'VALUE_ID',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub create {
    my $class = shift;
    my $bx = $class->define_boolexpr(@_);

    unless ($bx->value_for('value_class_name')) {
        my $pp_id = $bx->value_for('processing_profile_id');

        my $pp = Genome::ProcessingProfile->get($pp_id);
        die "invalid processing profile id $pp_id!" unless $pp;

        my $name = $bx->value_for('name');
        die "no name specified when constructing a processing profile param!" unless $name;

        my $pmeta = $pp->__meta__->property($name);
        die "no property $name found on processing profile " . $pp->__display_name__ unless $pmeta;

        my $value_class_name = $pmeta->_data_type_as_class_name;

        $bx = $bx->add_filter(value_class_name => $value_class_name);
    }

    return $class->SUPER::create($bx);
}

# this has the functionality of the old "value" accessor
# we wanted to ensure we were no longer dependent on it
# ..but the HTML view needs something generic which will work
sub _value_scalar_or_object {
    my $self = shift;
    my $name = $self->name;
    return $self->processing_profile->$name(@_);
}

sub value {
    Carp::confess("The system is calling value() on a Genome::ProcessingProfile::Param.  The old functionality of value() is not compatible with the new.  Code should go throuh the accessor on the processing profile, or call _value_scalar_or_object _IF_ it is internal to the profile")
}

1;

