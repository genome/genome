package Genome::Model::Tools::Analysis::Coverage;

use strict;
use warnings;

use Genome;

use XML::Simple;
use File::Basename;

class Genome::Model::Tools::Analysis::Coverage {
    is => ['Genome::Model::Tools::Analysis'],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools for analysis of sequence coverage (readcounts, etc)",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis coverage --help ...
EOS
}

sub help_detail {
    return <<EOS
EOS
}

1;
