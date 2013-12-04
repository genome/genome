package Genome::Config::Parser;

use warnings;
use strict;

class Genome::Config::Parser {
    is          => 'UR::Object',
    is_abstract => 1,
    doc         => 'abstract superclass for parsing various config formats',
};

sub parse {
    my $class = shift;
    my $file_path = shift;

    my $filetype = $class->_determine_filetype($file_path);
    my $concrete_parser = $class->_get_concrete_parser($filetype);

    return $concrete_parser->parse($file_path);
}

sub _determine_filetype {
    my $class = shift;
    my $file_path = shift;

    if ($file_path =~ /\.yml$/i || $file_path =~ /\.yaml$/i) {
        return 'yaml';
    } else {
        die("Unrecognized file type! $file_path");
    }
}

sub _get_concrete_parser {
    my $class = shift;
    my $filetype = shift;

    return sprintf('%s::%s', $class->class, uc($filetype));
}

1;
