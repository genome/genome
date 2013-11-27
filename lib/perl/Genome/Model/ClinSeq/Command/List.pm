package Genome::Model::ClinSeq::Command::List;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::ClinSeq::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::ClinSeq' 
        },
        show => { default_value => 'id,name,wgs_model,exome_model,tumor_rnaseq_model,normal_rnaseq_model' },
    ],
    doc => 'list clinseq genome models',
};

1;

