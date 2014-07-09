package Genome::Framework::Component::Reporter::SingleFile;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Framework::Component::Reporter::SingleFile {
    is => 'Genome::VariantReporting::Framework::Component::Reporter',
    has => [
        file_name => {},
    ],
    has_transient_optional => [
        _output_fh => {},
    ],
};

sub initialize {
    my $self = shift;
    $self->SUPER::initialize(@_);

    my $output_dir = shift;
    Genome::Sys->create_directory($output_dir);
    my $fh = Genome::Sys->open_file_for_writing(File::Spec->join($output_dir, $self->file_name));
    $self->_output_fh($fh);
    return;
}

sub finalize {
    my $self = shift;
    $self->SUPER::finalize(@_);

    $self->_output_fh->close;
    return;
}

1;
