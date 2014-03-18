package Genome::Model::GenotypeMicroarray::Command::List;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1, 
            value => 'Genome::Model::GenotypeMicroarray' 
        },
        show => { default_value => 'id,subject,dbsnp_build,genotype_vcf', },
    ],
    doc => 'List genotype microarray models',
};

1;

