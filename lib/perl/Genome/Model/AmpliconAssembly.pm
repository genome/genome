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

sub build_subclass_name {
    return 'amplicon-assembly';
}

1;

#$HeadURL$
#$Id$
