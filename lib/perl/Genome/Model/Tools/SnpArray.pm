package Genome::Model::Tools::SnpArray;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::SnpArray {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for working with array-based SNP genotype data.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt snp-array ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

