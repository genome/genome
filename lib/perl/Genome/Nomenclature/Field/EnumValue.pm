package Genome::Nomenclature::Field::EnumValue;

use strict;
use warnings;

use Command::Dispatch::Shell;
use Genome;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use JSON::XS;

class Genome::Nomenclature::Field::EnumValue {
    table_name => 'web.nomenclature_enum_value',
    id_by => [
        id => {
            is => 'Text',
            len => 255,
        },
    ],
    has => [
        nomenclature_field_id => {
            is => 'Text',
            len => 255,
            doc => 'Nomenclature field id',
        },
        nomenclature_field => {
            is => 'Genome::Nomenclature::Field',
            id_by => 'nomenclature_field_id',
            doc => 'Genome::Nomenclature field',
            constraint_name => 'GNEV_FK',
        },
        value => {
            is => 'Text',
            len => 255,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'Nomenclature enumerated field definitions',
};


1;
