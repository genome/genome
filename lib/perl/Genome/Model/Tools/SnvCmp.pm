package Genome::Model::Tools::SnvCmp;

use strict;
use warnings;

use Genome;

my $SNVCMP_DEFAULT = '1.0.0';

class Genome::Model::Tools::SnvCmp {
    is  => 'Command',
    is_abstract => 1,
    has_input => [
        use_version => {
            is  => 'Version', 
            doc => 'snvcmp version to be used.  default_value='. $SNVCMP_DEFAULT,
            is_optional   => 1, 
            default_value => $SNVCMP_DEFAULT,
        },
    ],
};


sub help_brief {
    "Tools to run snvcmp, a snv comparison tool.";
}

sub help_synopsis {
    "gmt snvcmp ...";
}

sub help_detail {                           
    "used to invoke snvcmp commands";
}

sub snvcmp_path {
    my $self = shift;
    return "/usr/bin/snvcmp";
}



1;

