package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath;

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Spec;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath { 
    is => 'Command::V2',
    is_abstract => 1,
    has_input => {
        source_path => {
            is => 'Text',
            doc => 'Source path of sequences to get.',
        },
       working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
    },
    has_output => {
        destination_path => {
            calculate_from => [qw/ working_directory source_path_basename /],
            calculate => q| return File::Spec->join($working_directory, $source_path_basename); |,
            doc => 'Final destination path.',
        }, 
        destination_md5_path => {
            calculate_from => [qw/ destination_path /],
            calculate => q| return Genome::InstrumentData::Command::Import::WorkFlow::Helpers->original_md5_path_for($destination_path); |,
            doc => 'Final destination MD5 path.',
        }, 
    },
    has_calculated => {
        source_path_basename => {
            calculate_from => [qw/ source_path /],
            calculate => q| return File::Basename::basename($source_path); |,
        },
    },
    has_constant_calculated => {
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    },
};

sub execute {
    my $self = shift;

    my $retrieve_source_path = $self->retrieve_source_path;
    return if not $retrieve_source_path;

    my $retrieve_source_md5 = $self->retrieve_source_md5;
    return if not $retrieve_source_md5;

    return 1;
}

sub retrieve_source_path {
    my $self = shift;
    $self->debug_message('Retrieve source path...');

    my $source_path = $self->source_path;
    $self->debug_message("Source path: $source_path");
    my $source_path_sz = $self->_source_path_size;
    $self->debug_message("Source path size: ".(defined $source_path_sz ? $source_path_sz : 'NA'));
    if ( not defined $source_path_sz or $source_path_sz == 0 ) { # error that the file does not exist
        $self->error_message('Source file does not have any size!');
        return;
    }

    my $destination_path = $self->destination_path;
    $self->debug_message("Destination path: $destination_path");

    my $retrieve_ok = $self->_retrieve_source_path;
    return if not $retrieve_ok;

    my $destination_path_sz = -s $destination_path;
    $self->debug_message("Destination path size: $destination_path_sz");
    if ( not $source_path_sz or $source_path_sz != $destination_path_sz ) {
        $self->error_message("Retrieve succeeded, but source path size is different from destination path! $source_path_sz <=> $destination_path_sz");
        return;
    }

    $self->debug_message('Retrieve source path...done');
    return 1;
}

sub retrieve_source_md5 {
    my $self = shift;
    $self->debug_message('Create destination MD5 path...');

    my $source_md5 = $self->source_md5;
    if ( not $source_md5) {
        $self->debug_message('Source MD5 is not available! It will not be saved.');
        return 1;
    }

    $self->debug_message("Source MD5: $source_md5");
    if ( $source_md5 !~ /^[0-9a-f]{32}$/i ) {
        $self->debug_message('Invalid source MD5! It will not be saved.');
        return 1;
    }

    my $destination_md5_path = $self->destination_md5_path;
    $self->debug_message("Destination MD5 path: $destination_md5_path");

    my $fh = Genome::Sys->open_file_for_writing($destination_md5_path);
    $fh->say( join('  ', $source_md5, $self->destination_path) );
    $fh->close;

    $self->debug_message('Create destination MD5 path...done');
    return 1;
}

1;

