package Genome::Model::Tools::Solexa;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;

class Genome::Model::Tools::Solexa {
    is => ['Command'],
    has_optional => [
                     version => {
                                 is    => 'string',
                                 doc   => 'version of Solexa application to use',
                             },
                     _tmp_dir => {
                                  is => 'string',
                                  doc => 'a temporary directory for storing files',
                              },
                 ]
};

sub help_brief {
    "tools to work with 454 reads"
}

sub help_detail {                           # This is what the user will see with --help <---
    return <<EOS

EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    unless (Genome::Config->arch_os =~ /64/) {
        $self->error_message('Most Solexa tools must be run from 64-bit architecture');
        return;
    }
    my $tempdir = File::Temp::tempdir(CLEANUP => 1);
    $self->_tmp_dir($tempdir);

    return $self;
}

sub resolve_454_path {
    return $ENV{GENOME_SW} . '/454/';
}

sub bin_path {
    my $self = shift;
    return $self->resolve_454_path .'offInstrumentApps-'. $self->version .'/bin';
}

1;

