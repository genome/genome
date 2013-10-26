package Genome::SoftwareResult::Input;

use strict;
use warnings;

use Genome;
class Genome::SoftwareResult::Input {
    table_name => 'result.input',
    type_name => 'software result input',
    id_by => [
        software_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'software_result_id',
        },
        name => {
            is => 'Text',
            len => 255,
            column_name => 'input_name',
        },
    ],
    has => [
        value_class_name => {
            is => 'Text',
            len => 255,
            is_optional => 1,
        },
        value_id => {
            is => 'Text',
            len => 1000,
            column_name => 'input_value',
        },
        value_obj => {
            is => 'UR::Object',
            id_by => 'value_id',
            id_class_by => 'value_class_name',
        },
        # after the new API is released and old snapshots go away, invert the column assingments with those above
        _new_name => {
            is => 'Text',
            len => 255,
            column_name => 'name',
            is_optional => 1,
        },
        _new_value => {
            is => 'Text',
            len => 1000,
            column_name => 'value_id',
            is_optional => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub __display_name__ {
    my $self = shift;
    my $name = $self->name;
    my $sr = $self->software_result;
    my @value;
    my $value_class_name = $self->value_class_name;
    my $value_id = $self->value_id;
    my $value;
    if ($value_class_name) {
        $value = $value_class_name->get($value_id);
    }
    else {
        $value = UR::Value::Text->get($value_id);
    }
    return "$name:" . $value->__display_name__;
}

# this will sync up new columns with the old ones
# once all old snapshots are gone, we will switch to the new columns and remove this

sub create {
    my $class = shift;
    my $bx = $class->define_boolexpr(@_);

    unless ($bx->value_for('value_class_name')) {
        my $sr_id = $bx->value_for('software_result_id');

        my $sr = Genome::SoftwareResult->get($sr_id);
        die "invalid software result id $sr_id!" unless $sr;

        my $name = $bx->value_for('name');
        $name =~ s/-.*//;
        die "no name specified when constructing a software result input!" unless $name;

        my $pmeta = $sr->__meta__->property($name);
        die "no property $name found on software result " . $sr->__display_name__ unless $pmeta;

        my $value_class_name = $pmeta->_data_type_as_class_name;

        $bx = $bx->add_filter(value_class_name => $value_class_name);
    }

    my $self = $class->SUPER::create($bx);
    return unless $self;
    $self->_new_name($self->name);
    $self->_new_value($self->value_id);
    return $self;
}

# this has the functionality of the old "value" accessor
# we wanted to ensure we were no longer dependent on it
# ..but the HTML view needs something generic which will work
sub _value_scalar_or_object {
    my $self = shift;
    my $name = $self->name;
    return $self->software_result->$name(@_);
}

sub value {
    Carp::confess("The system is calling value() on a Genome::ProcessingProfile::input.  The old functionality of value() is not compatible with the new.  Code should go throuh the accessor on the processing profile, or call _value_scalar_or_object _IF_ it is internal to the profile")
}

sub input_name {
    Carp::confess("using input_name!");
}

sub input_value {
    Carp::confess("using input_value!");
}

1;

