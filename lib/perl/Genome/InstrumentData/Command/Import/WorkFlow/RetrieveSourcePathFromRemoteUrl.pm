package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl;

use strict;
use warnings;

use Genome;

use LWP::UserAgent;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromRemoteUrl { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
};

sub _retrieve_path {
    my ($self, $from, $to) = @_;
    $self->debug_message('Retrieve source path via http...');

    $self->debug_message("From: $from");
    $self->debug_message("To: $to");

    my $agent = LWP::UserAgent->new;
    my $response = $agent->get($from, ':content_file' => $to);
    if ( not $response->is_success ) {
        $self->error_message($response->message) if $response->message;
        $self->error_message('GET failed for remote path!');
        return
    }

    my $from_sz = $response->headers->content_length;
    $self->debug_message("From size: $from_sz");
    my $to_sz = -s $to;
    $self->debug_message("To size: $to_sz");
    if ( not $to_sz or $to_sz != $from_sz ) {
        $self->error_message("GET remote path succeeded, but destination size is different from original! $from_sz vs $to_sz");
        return;
    }

    $self->debug_message('Retrieve source path from remote url...done');
    return 1;
}
1;

