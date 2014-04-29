package Genome::Annotation::Reporter::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Reporter::Base {
    is => 'Genome::Annotation::ComponentBase',
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
    my $fh = Genome::Sys->open_file_for_writing(File::Spec->join($output_dir, $self->file_name));
    $self->_output_fh($fh);
}

sub finalize {
    my $self = shift;
    $self->_output_fh->close;
}

1;
