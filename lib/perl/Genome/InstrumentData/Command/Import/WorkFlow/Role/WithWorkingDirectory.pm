package Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory;

use strict;
use warnings;

require File::Basename;
require File::Spec;

use Genome;

role Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory {
    has => [
        working_directory => {
            is => 'DirectoryPath',
            is_optional => 0,
            is_input => 1,
            doc => 'Working directory to put output files.',
        },
    ],
};

sub get_working_path_for_file_path {
    return File::Spec->join( $_[0]->working_directory, File::Basename::basename($_[1]) );
}

sub get_working_bam_path_with_new_extension {
    my ($self, $bam_path, @ext) = @_;

    if ( not $bam_path =~ s/\.bam$// ) {
        $self->fatal_message('Cannot insert new extension into bam path because it does not end with .bam! %', $bam_path);
    }

    return join('.', $self->get_working_path_for_file_path($bam_path), @ext, 'bam');
}

1;

