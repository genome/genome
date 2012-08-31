package Genome::Model::Tools::Allpaths::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Allpaths::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
	    version => {
            is => 'Text',
            doc => 'Version of ALLPATHS to use',
            valid_values => [ '41055' ],
            default_value => "41055",
        },
    ],
};

sub executable_for_version {
    my ($self, $name) = @_;
    return $name."-".$self->version;
}

1;

