package Genome::Cromwell;

use strict;
use warnings;

use JSON qw(to_json from_json);
use HTTP::Request;
use HTTP::Request::Common qw(POST);
use LWP::UserAgent;
use IO::Socket::SSL qw();
use IPC::Run qw();
use Scope::Guard;

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
            is_calculated => 1,
            calculate => q{ Genome::Config::get('cromwell_api_server') },
        },
        user_agent => {
            is => 'LWP::UserAgent',
            is_constant => 1,
            calculate => q{
                my $ua = LWP::UserAgent->new();
                $ua->ssl_opts( verify_hostname => 0,  SSL_verify_mode => IO::Socket::SSL::SSL_VERIFY_NONE );
                return $ua;
            },
        },
        ],
};

sub status {
    my $class = shift;
    my $workflow_id = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url($workflow_id, 'status');

    return $self->_make_json_request('GET', $url);
}

sub outputs {
    my $class = shift;
    my $workflow_id = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url($workflow_id, 'outputs');

    return $self->_make_json_request('GET', $url);
}

sub metadata {
    my $class = shift;
    my $workflow_id = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url($workflow_id, 'metadata');
    $url .= '?excludeKey=submittedFiles&excludeKey=inputs&excludeKey=outputs&expandSubWorkflows=false';

    return $self->_make_json_request('GET', $url);
}

sub query {
   my $class = shift;
   my $query_options = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url('query');

    return $self->_make_json_request('POST', $url, $query_options);
}

sub timing {
    my $class = shift;
    my $workflow_id = shift;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url($workflow_id, 'timing');

    my $content = $self->_make_request('GET', $url);
    return $content;
}

sub submit_workflow {
    my $class = shift;
    my ($definition, $inputs, $dependencies, $options) = @_;

    my $self = $class->_singleton_object;
    my $url = $self->_request_url();
    my $content = [ workflowSource  => [$definition],
                    workflowInputs => [$inputs],
                    workflowDependencies => [$dependencies],
                    workflowOptions => [$options] ];
    my $req = POST( $url,
                    Content_Type => 'multipart/form-data',
                    Content => $content );

    return from_json($self->_send_request($req));
}


sub _request_url {
    my $class = shift;
    my @parts = @_;

    my $self = $class->_singleton_object;

    my $url = $self->server_url;
    $url .= join('/', 'api', 'workflows', $self->api_version, @parts);

    return $url;
}

sub _make_json_request {
    my $class = shift;

    my $content = $class->_make_request(@_);
    return from_json($content);
}

sub _send_request {
    my $class = shift;
    my $req = shift;

    my $self = $class->_singleton_object;
    ATTEMPT: for (1..5) {
        my $response = $self->user_agent->request($req);

        unless ($response->is_success) {
            $self->fatal_message('Error querying server: %s', $response->status_line);
            sleep 5;
            next ATTEMPT;
        }

        return $response->decoded_content;
    }

    $self->fatal_message('Failed to query server after serveral attempts.');
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
    return $self->_send_request(HTTP::Request->new(@request_args));
}

sub cromwell_jar_cmdline {
    my $class = shift;
    my $config = shift;

    my $truststore_file = Genome::Config::get('cromwell_truststore_file');
    my $truststore_auth = Genome::Config::get('cromwell_truststore_auth');

    my @cmd = (
        '/usr/bin/java',
        sprintf('-Dconfig.file=%s', $config),
        sprintf('-Djavax.net.ssl.trustStorePassword=%s', $truststore_auth),
        sprintf('-Djavax.net.ssl.trustStore=%s', $truststore_file),
        '-jar', '/opt/cromwell.jar',
    );
    return @cmd;
}

sub spawn_local_server {
    my $class = shift;
    my $config = shift;

    my @jar_cmdline = $class->cromwell_jar_cmdline($config);

    $class->debug_message('Spawning local cromwell server: %s', join(' ', @jar_cmdline));

    my $in = '';
    my $out = '';

    my $harness = IPC::Run::start([@jar_cmdline, 'server'], \$in, \$out);
    my $env_guard = Genome::Config::set_env('cromwell_api_server', 'http://localhost:8000/');

    $harness->pump until $out =~ /Cromwell [^ ]+ service started on/;

    my $guard_closure = sub {
        $env_guard = undef;
        $harness->kill_kill();
    };

    my $server_guard = Scope::Guard->new($guard_closure);
    return $server_guard;
}


1;
