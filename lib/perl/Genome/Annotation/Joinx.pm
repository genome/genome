package Genome::Annotation::Joinx;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Joinx {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
        variant_type => { 
            is => 'Text', 
        },
        info_string => {
            is => 'Text',
        },
        version => {
            is => 'Text',
        },
    ],
    has_optional_output => [
        software_result => {
            is => 'Genome::Annotation::Joinx::Result',
            doc => 'The software result created during command execution',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->software_result(Genome::Annotation::Joinx::Result->get_or_create($self->input_hash));
    return 1;
}
