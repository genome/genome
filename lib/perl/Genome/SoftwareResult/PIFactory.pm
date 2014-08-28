package Genome::SoftwareResult::PIFactory;

use strict;
use warnings;

use Genome;

use Carp qw(croak);

sub import {
    my $classname = shift;

    my $subclassname = (caller)[0];
    my $type = lc((split('::', $subclassname))[-1]);

    my $subclass = UR::Object::Type->define(
        class_name => $subclassname,
        isa => 'Genome::SoftwareResult::PIBase',
        table_name => 'result.' . $type,
        id_by => [
            software_result_id => {
                is => 'Text',
                column_name => 'software_result_id',
            },
            name => {
                is => 'Text',
                len => 255,
                column_name => $type . '_name',
            },
        ],
        has => [
            value_id => {
                is => 'Text',
                len => 1000,
                column_name => $type . '_value',
            },
            value_class_name => {
                is => 'Text',
                len => 255,
                is_optional => 1,
            },
            value_obj => {
                is => 'UR::Object',
                id_by => 'value_id',
                id_class_by => 'value_class_name',
            },
            software_result => {
                is => 'Genome::SoftwareResult',
                id_by => 'software_result_id',
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
    );
}

1;
