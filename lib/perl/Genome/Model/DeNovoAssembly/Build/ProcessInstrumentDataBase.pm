package Genome::Model::DeNovoAssembly::Build::ProcessInstrumentDataBase;

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Build::ProcessInstrumentDataBase {
    is => 'Command::V2',
    has_input => {
        build => {
            is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1,
        },
    },
    has_optional_constant => {
        lsf_resource => {
            default_value => "-R 'select[mem>32000 && gtmp>200] rusage[mem=32000:gtmp=200] span[hosts=1]' -M 32000000",
        },
    },
};

sub lsf_resource { # in 2 places b/c workflow requires a 'default_value' for this property
    $_[0]->__meta__->property_meta_for_name('lsf_resource')->default_value;
}

1;

