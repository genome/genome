package Genome::InstrumentDataAttribute;

use strict;
use warnings;
use Genome;

class Genome::InstrumentDataAttribute {
    table_name => 'INSTRUMENT_DATA_ATTRIBUTE',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Represents a particular attribute of an instrument data object',
    id_by => [
        instrument_data_id => {
            is => 'Text',
        },
        attribute_label => {
            is => 'Text',
        },
        attribute_value => {
            is => 'Text',
        },
        nomenclature => {
            is => 'Text',
        },
    ],
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id',
        },
    ],
};

sub create {
    my $class = shift;
    my $bx = $class->define_boolexpr(@_);
    # TODO This is a workaround that allows nomenclature to be in the id_by block
    # and have a default value. Doing so in the class definition doesn't work due
    # to some sort of UR bug that Tony is aware of.
    unless ($bx->specifies_value_for('nomenclature')) {
        $bx = $bx->add_filter('nomenclature' => 'WUGC');
    }
    return $class->SUPER::create($bx);
}
1;

