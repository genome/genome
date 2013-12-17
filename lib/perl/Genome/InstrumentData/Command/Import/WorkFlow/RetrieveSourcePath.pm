package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath;

use strict;
use warnings;

use Genome;

require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath { 
    is => 'Command::V2',
    has_input => [
        source_path => {
            is => 'Text',
            doc => 'Source path of sequences to get.',
        },
       working_directory => {
            is => 'Text',
            doc => 'Detination directory for source path.',
        },
    ],
    has_output => [
        destination_path => {
            calculate_from => [qw/ working_directory source_path /],
            calculate => q| return $self->working_directory.'/'.File::Basename::basename($source_path); |,
            doc => 'Final destination path.',
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

    my $retrieve_source_path = $self->_retrieve_source_path;
    return if not $retrieve_source_path;

    my $retrieve_source_md5 = $self->_retrieve_source_md5;
    return if not $retrieve_source_md5;

    return 1;
}

sub _retrieve_source_path {
    my $self = shift;
    $self->status_message('Retrieve source path...');

    my $retrieve_ok = $self->retrieve_path($self->source_path, $self->destination_path);
    return if not $retrieve_ok;
    
    $self->status_message('Retrieve source path...done');
    return 1;
}

sub _retrieve_source_md5 {
    my $self = shift;
    $self->status_message('Retrieve source MD5 path...');

    my $md5_path = $self->source_path.'.md5';
    my $md5_size = $self->helpers->file_size($md5_path);
    if ( not $md5_size ) {
        $self->status_message('Source MD5 path does not exist...skip');
        return 1;
    }

    my $retrieve_ok = $self->retrieve_path($md5_path, $self->destination_path.'.md5-orig');
    return if not $retrieve_ok;

    $self->status_message('Retrieve source MD5 path...done');
    return 1;
}

sub retrieve_path {
    my ($self, $from, $to) = @_;

    Carp::confess('No from path to retrieve file!') if not $from;
    Carp::confess('No to path to retrieve file!') if not $to;

    if ( $from =~ /^http/ ) {
        return $self->_retrieve_remote_path($from, $to);
    }
    else {
        return $self->_retrieve_local_path($from, $to);
    }
}

sub _retrieve_local_path {
    my ($self, $from, $to) = @_;
    $self->status_message('Retrieve local path...');

    $self->status_message("From: $from");
    $self->status_message("To: $to");
    my $move_ok = File::Copy::copy($from, $to);
    if ( not $move_ok ) {
        $self->error_message('Copy failed!');
        return;
    }

    my $from_sz = -s $from;
    $self->status_message("From size: $from_sz");
    my $to_sz = -s $to;
    $self->status_message("To size: $to_sz");
    if ( not $to_sz or $to_sz != $from_sz ) {
        $self->error_message("Copy succeeded, but destination size is diffeerent from original! $to_sz vs $from_sz");
        return;
    }

    $self->status_message('Retrieve local path...done');
    return 1;
}

sub _retrieve_remote_path {
    my ($self, $from, $to) = @_;
    $self->status_message('Retrieve remote path...');

    $self->status_message("From: $from");
    $self->status_message("To: $to");

    my $agent = LWP::UserAgent->new;
    my $response = $agent->get($from, ':content_file' => $to);
    if ( not $response->is_success ) {
        $self->error_message($response->message) if $response->message;
        $self->error_message('GET failed for remote path!');
        return
    }

    my $from_sz = $response->headers->content_length;
    $self->status_message("From size: $from_sz");
    my $to_sz = -s $to;
    $self->status_message("To size: $to_sz");
    if ( not $to_sz or $to_sz != $from_sz ) {
        $self->error_message("GET remote path succeeded, but destination size is different from original! $to_sz vs $from_sz");
        return;
    }

    $self->status_message('Retrieve remote path...done');
    return 1;
}
1;

