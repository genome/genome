package Genome::InstrumentDataAttribute;

use strict;
use warnings;
use Genome;

class Genome::InstrumentDataAttribute {
    table_name => 'instrument.data_attribute',
    id_by => [
        instrument_data_id => {
            is => 'Text',
            len => 64,
        },
        attribute_label => {
            is => 'Text',
            len => 64,
        },
        attribute_value => {
            is => 'Text',
            len => 512,
        },
        nomenclature => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id',
            constraint_name => 'IDA_ID_FK',
        },
        # TODO: we have been simplifying the name/value stuff for some time
        # Switch to these by default, and test the inversion.
        name => {
            via => '__self__',
            to => 'attribute_label',
        },
        value => {
            via => '__self__',
            to => 'attribute_value',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Represents a particular attribute of an instrument data object',
};

sub create {
    my $class = shift;
    my $bx = $class->define_boolexpr(@_);
    # TODO This is a workaround that allows nomenclature to be in the id_by block
    # and have a default value. Doing so in the class definition doesn't work due
    # to some sort of UR bug that Tony is aware of.
    unless ($bx->specifies_value_for('nomenclature')) {
        $bx = $bx->add_filter('nomenclature' => $ENV{GENOME_NOMENCLATURE_DEFAULT});
    }
    return $class->SUPER::create($bx);
}

1;

