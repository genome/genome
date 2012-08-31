package Genome::Nomenclature::Field::EnumValue;

use strict;
use warnings;

use Command::Dispatch::Shell;
use Genome;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use JSON::XS;

class Genome::Nomenclature::Field::EnumValue {
    table_name => 'GENOME_NOMENCLATURE_ENUM_VALUE',
    id_generator => '-uuid',
    id_by => {
        'id' => {is=>'Text', len=>64}
    },
    has => [
        nomenclature_field_id => {
            is=>'Text', 
            len=>255, 
            doc => 'Nomenclature field id'
        },
        nomenclature_field => {
            is=>'Genome::Nomenclature::Field', 
            doc => 'Genome::Nomenclature field',
            id_by => 'nomenclature_field_id'
        },
        value => {
            is=>'Text',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Nomenclature enumerated field definitions'
};


1;
