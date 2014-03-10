package Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs;

use strict;
use warnings;

use Genome;

use Archive::Extract;
$Archive::Extract::PREFER_BIN = 1;
require File::Basename;
require File::Find;
require File::Path;

class Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs { 
    is => 'Command::V2',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'Destination directory for fastqs.',
        },
        archive_path => { 
            is => 'Text',
            doc => 'Path of the fastq archive to unpack.',
        },
    ],
    has_output => [ 
        fastq_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'The paths of the extracted fastqs.',
        },
    ],
    has_optional => [
        extract_directory => {
            calculate_from => [qw/ working_directory /],
            calculate => q( return $working_directory.'/extract'; ),
        },
        extracted_fastqs => {
            is => 'Array',
        },
    ],
    has_constant_calculated => [
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Archive to Fastqs...');

    my $archive_to_fastqs_ok = $self->_archive_to_fastqs;
    return if not $archive_to_fastqs_ok;

    my $move_fastqs = $self->_move_extracted_fastqs_to_working_directory;
    return if not $move_fastqs;

    my $cleanup_ok = $self->_cleanup;
    return if not $cleanup_ok;

    $self->debug_message('Archive to Fastqs...done');
    return 1;
}

sub _archive_to_fastqs {
    my $self = shift;
    $self->debug_message('Extract...');

    my $archive_path = $self->archive_path;
    $self->debug_message("Archive path: $archive_path");

    my $extractor = Archive::Extract->new(archive => $self->archive_path);
    if ( not $extractor ) {
        $self->error_message('Failed to create extractor!');
        return;
    }

    my $extract_directory = $self->extract_directory;
    mkdir $extract_directory;

    my $extract_ok = $extractor->extract(to => $extract_directory);
    if ( not $extract_ok ) {
        my $error_message = $extractor->error;
        $self->error_message($error_message) if defined $error_message;
        $self->error_message('Failed to extract!');
        return;
    }

    my @extracted_fastqs = grep { /\.f(ast)?q$/ } @{$extractor->files};
    if ( not @extracted_fastqs ) {
        $self->error_message('No fastqs found in archive!');
        return;
    }
    $self->extracted_fastqs(\@extracted_fastqs);


    $self->debug_message('Extract...done');
    return 1;
}

sub _move_extracted_fastqs_to_working_directory {
    my $self = shift;
    $self->debug_message('Move extracted fastqs to working directory...');

    my @fastq_paths;
    for my $extracted_fastq ( @{$self->extracted_fastqs} ) {
        my $extracted_fastq_path = $self->extract_directory.'/'.$extracted_fastq;
        $self->debug_message("Extracted FASTQ path: ".$extracted_fastq_path);
        my $extracted_fastq_basename = File::Basename::basename($extracted_fastq_path);
        my $fastq_path = $self->working_directory.'/'.$extracted_fastq_basename;
        my $move_ok = $self->helpers->move_path($extracted_fastq_path, $fastq_path);
        Carp::confess('Failed to move file!') if not $move_ok;
        $self->debug_message("FASTQ path: ".$fastq_path);
        push @fastq_paths, $fastq_path;
    }

    if ( not @fastq_paths ) {
        $self->error_message('No fastqs found in archive!');
        return;
    }
    $self->fastq_paths(\@fastq_paths);

    $self->debug_message('Move extracted fastqs to working directory...done');
    return 1;
}

sub _cleanup {
    my $self = shift;
    $self->helpers->remove_paths_and_auxiliary_files($self->archive_path);
    File::Path::remove_tree($self->extract_directory);
    return 1;
}

1;

