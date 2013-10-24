package Genome::InstrumentData::Command::Import::WorkFlow::ArchiveToFastqs;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
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
            doc => 'The paths of the unarchived fastqs.',
        },
    ],
    has_optional_ => [
        unarchive_directory => {
            calculate_from => [qw/ working_directory /],
            calculate => q( return $working_directory.'/unarchive'; ),
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
    $self->status_message('Archive to Fastqs...');

    my $archive_to_fastqs_ok = $self->_archive_to_fastqs;
    return if not $archive_to_fastqs_ok;

    my $find_fastqs = $self->_find_fastqs;
    return if not $find_fastqs;

    my $cleanup_ok = $self->_cleanup;
    return if not $cleanup_ok;

    $self->status_message('Archive to Fastqs...done');
    return 1;
}

sub _archive_to_fastqs {
    my $self = shift;
    $self->status_message('Unarchive...');

    my $archive_path = $self->archive_path;
    $self->status_message("Archive path: $archive_path");
    my $unarchive_directory = $self->unarchive_directory;
    mkdir $unarchive_directory;
    my $cmd = "tar xzf $archive_path --directory=$unarchive_directory";

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run tar!');
        return;
    }

    $self->status_message('Unarchive...done');
    return 1;
}

sub _find_fastqs {
    my $self = shift;
    $self->status_message('Find fastqs...');
    
    my @fastq_paths;
    File::Find::find(
        sub{
            return if not /\.fa?s?t?q$/;
            $self->status_message("Unarchived FASTQ path: ".$File::Find::name);
            my $fastq_path = $self->working_directory.'/'.$_;
            my $move_ok = $self->helpers->move_path($File::Find::name, $fastq_path);
            Carp::confess('Failed to move file!') if not $move_ok;
            $self->status_message("FASTQ path: ".$fastq_path);
            push @fastq_paths, $fastq_path;
        },
        $self->working_directory,
    );
    $self->fastq_paths(\@fastq_paths);

    $self->status_message('Find fastqs...done');
    return 1;
}

sub _cleanup {
    my $self = shift;
    $self->helpers->remove_paths_and_auxiliary_files($self->archive_path);
    File::Path::remove_tree($self->unarchive_directory);
    return 1;
}

1;

