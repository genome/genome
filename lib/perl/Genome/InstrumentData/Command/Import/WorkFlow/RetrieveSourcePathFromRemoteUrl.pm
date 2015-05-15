package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl;

use strict;
use warnings;

use Genome;

use LWP::UserAgent;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};

sub _retrieve_source_path {
    my $self = shift;

    my $agent = LWP::UserAgent->new;
    my $response = $agent->get($self->source_path, ':content_file' => $self->destination_path);
    if ( not $response->is_success ) {
        $self->error_message($response->message) if $response->message;
        $self->error_message('Remote GET from %s to %s to failed!', $self->source_path, $self->destination_path);
        return;
    }

    return 1;
}

sub _load_source_md5 {
    my $self = shift;

    return if not $self->source_file->md5_path_size;
    my $agent = LWP::UserAgent->new;
    my $response = $agent->get($self->source_file->md5_path);
    if ( not $response->is_success ) {
        return;
    }

    my ($md5) = split(/\s+/, $response->decoded_content);
    return $md5;
}

1;

