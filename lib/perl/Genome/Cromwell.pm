package Genome::Cromwell;

use strict;
use warnings;

use JSON qw(to_json from_json);
use HTTP::Request;
use LWP::UserAgent;

use Genome;

class Genome::Cromwell {
    is => 'UR::Singleton',
    has => [
        api_version => {
            is => 'Text',
            is_constant => 1,
            value => 'v1',
        },
        server_url => {
            is => 'Text',
            is_constant => 1,
            value => Genome::Config::get('cromwell_api_server'),
        },
        user_agent => {
            is => 'LWP::UserAgent',
            is_constant => 1,
            calculate => q{ LWP::UserAgent->new(); },
        },
    ],
};


sub outputs {
    my $class = shift;
    my $workflow_id = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url($workflow_id, 'outputs');

    return $self->_make_request('GET', $url);
}

sub query {
   my $class = shift;
   my $query_options = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url('query');

    return $self->_make_request('POST', $url, $query_options);
}

sub _request_url {
    my $class = shift;
    my @parts = @_;

    my $self = $class->_singleton_object;

    my $url = $self->server_url;
    $url .= join('/', 'api', 'workflows', $self->api_version, @parts);

    return $url;
}

sub _make_request {
    my $class = shift;
    my ($method, $url, $data) = @_;

    my $self = $class->_singleton_object;
    my @headers = (
        'Accept' => 'application/json',
        'Accept-Encoding' => 'gzip',
    );

    my @request_args = ($method, $url);
    if (defined $data) {
        my $content = to_json($data);
        push @headers, 'Content-Type' => 'application/json';
        push @headers, 'Content-Length' => length $content;
        push @request_args, \@headers, $content;
    } else {
        push @request_args, \@headers;
    }

    my $req = HTTP::Request->new(@request_args);
    my $response = $self->user_agent->request($req);

    unless ($response->is_success) {
        $self->fatal_message('Error querying server: %s', $response->status_line);
    }

    my $content = $response->decoded_content;
    return from_json($content);
}

1;
