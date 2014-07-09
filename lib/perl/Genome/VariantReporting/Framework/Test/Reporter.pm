package Genome::VariantReporting::Framework::Test::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::Reporter {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    has => [
        file_name => {},
    ],
    has_transient_optional => [
        _output_fh => {},
    ],
};

sub name {
    return '__test__';
}

sub requires_interpreters {
    return qw(__test__);
}

sub initialize {
    my $self = shift;
    my $output_dir = shift;
    Genome::Sys->create_directory($output_dir);
    my $fh = Genome::Sys->open_file_for_writing(File::Spec->join($output_dir, $self->file_name));
    $self->_output_fh($fh);
}

sub finalize {
    my $self = shift;
    $self->_output_fh->close;
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    $self->_output_fh->print(pp($interpretations) . "\n");
}


1;
