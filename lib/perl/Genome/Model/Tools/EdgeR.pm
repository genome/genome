package Genome::Model::Tools::EdgeR;

use strict;
use warnings;

use Genome;
use Carp qw/confess/;
use Sys::Hostname;

my $DEFAULT_VER = '1.6';
my $MINIMUM_VER_FOR_RLIB = 1.5;

class Genome::Model::Tools::EdgeR {
    is  => 'Command',
    is_abstract => 1,
};
1;

