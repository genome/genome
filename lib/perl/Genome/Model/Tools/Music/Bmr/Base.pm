package Genome::Model::Tools::Music::Bmr::Base;

use strict;
use warnings;
use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Bmr::Base {
    is => ['Genome::Model::Tools::Music::Base'],
    is_abstract => 1,
    attributes_have => [
        file_format => {
            is => 'Text',
            is_optional => 1,
        }
    ],
    doc => "Mutational Significance In Cancer (BMR Calculations)"
};

sub _doc_authors {
    return " Cyriac Kandoth, Ph.D.";
}

sub _doc_see_also {
    return <<EOS
B<genome-music-bmr>(1),
B<genome-music>(1),
B<genome>(1)
EOS
}

1;
