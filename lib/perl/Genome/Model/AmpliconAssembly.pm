package Genome::Model::AmpliconAssembly;

use strict;
use warnings;

use Genome;

class Genome::Model::AmpliconAssembly {
    is => 'Genome::ModelDeprecated',
    has => [
    map({
            $_ => {
                via => 'processing_profile',
            }
        } Genome::ProcessingProfile::AmpliconAssembly->params_for_class
    ),
    ],
};

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

sub is_deprecated { 1 }

sub build_subclass_name {
    return 'amplicon-assembly';
}

1;
