package Genome::Model::Tools::Beagle;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;

class Genome::Model::Tools::Beagle {
    is => ['Command'],
    has_optional => [
                     version => {
                                 is    => 'string',
                                 doc   => 'version of Beagle application to use',
                             },
                     _tmp_dir => {
                                  is => 'string',
                                  doc => 'a temporary directory for storing files',
                              },
                 ]
};

sub help_brief {
    "Tools to run Beagle genetic analysis"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    my $tempdir = File::Temp::tempdir(CLEANUP => 1);
    $self->_tmp_dir($tempdir);

    return $self;
}

sub path_to_binary {
    return("java -Xmx14000m -jar $ENV{GENOME_SW}/beagle/installed/beagle.jar");
}
1;

