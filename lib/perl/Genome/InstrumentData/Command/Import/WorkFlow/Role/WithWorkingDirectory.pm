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

sub get_working_bam_path_with_new_extension {
    my ($self, $bam_path, @ext) = @_;

    my $bam_basename = File::Basename::basename($bam_path);
    if ( not $bam_basename =~ s/\.bam$// ) {
        $self->fatal_message('Cannot insert new extension into bam path because it does not end with .bam! %', $bam_path);
    }

    return File::Spec->join( $self->working_directory, join('.', $bam_basename, @ext, 'bam') );
}

1;

