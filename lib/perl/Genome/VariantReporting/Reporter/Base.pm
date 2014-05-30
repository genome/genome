package Genome::VariantReporting::Reporter::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Reporter::Base {
    is => 'Genome::VariantReporting::ComponentBase',
    is_abstract => 1,
    has => [
        file_name => {
            is => 'Text',
        },
    ],
    has_transient_optional => [
        _output_fh => {},
    ],
};

sub name {
    die "abstract";
}

sub requires_interpreters {
    die "abstract - must return a list of one or more interpreter names";
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

1;
