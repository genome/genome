package Genome::SoftwareResult::ImportedFile;

use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::ImportedFile {
    is => 'Genome::SoftwareResult::StageableSimple::SingleFile',
    has_input => [
        source_file_path => {
            is => 'Path',
        }
    ],
};

sub _run{
    my $self = shift;
    Genome::Sys->copy_file($self->source_file_path, $self->_temp_staging_file_path);
}

sub _file_name {
    my $self = shift;
    return File::Basename::basename($self->source_file_path);
}

1;
