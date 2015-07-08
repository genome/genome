package Genome::Config::Parser::YAML;

use warnings;
use strict;

use YAML::Syck;
use Params::Validate qw(validate_pos);

class Genome::Config::Parser::YAML {
    is => 'Genome::Config::Parser',
    doc => 'handles the parsing of a YAML config file'
};

sub parse {
    my ($self, $filename) = validate_pos(
        @_,
        1,
        {
            callbacks => {
                'File exists' => sub { -e $_[0]; },
                'File has size' => sub { -s $_[0]; },
                'File is YAML' => sub { $_[0] =~ /.+\.(yaml|yml)$/i },
            }
        }
    );
    return LoadFile($filename);
}

sub write {
    my ($self, $filename, $data) = validate_pos(@_, 1, 1, 1);
    return DumpFile($filename, $data);
}

1;
