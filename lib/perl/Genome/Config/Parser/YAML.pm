package Genome::Config::Parser::YAML;

use warnings;
use strict;

use YAML::Syck;

class Genome::Config::Parser::YAML {
    is => 'Genome::Config::Parser',
    doc => 'handles the parsing of a YAML config file'
};

sub parse {
    my ($self, $filename) = @_;
    die('Failed to provide filename!') unless $filename;
    die("$filename doesn't exist!") unless (-e $filename);
    die("$filename is empty!") unless (-s $filename);
    die("$filename doesn't appear to be a YAML file!") unless $filename =~ /.+\.(yaml|yml)$/i;
    return LoadFile($filename);
}

1;