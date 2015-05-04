package Genome::Model::Tools::CgHub::Query;

use strict;
use warnings;

use Genome;

use File::Temp;
use HTTP::Request;
use LWP::UserAgent;

class Genome::Model::Tools::CgHub::Query {
    is => "Genome::Model::Tools::CgHub::Base",
    has_input => {
        query => {
            is => "Text",
            doc => 'The query to send to CG Hub.'
        },
    },
    has_optional_output => {
        xml_file => {
            is => "Text",
            doc => 'Save metadata XML output to this file.',
        },
    },
    has_optional_transient_output => {
        metadata_xml => { is =>'Text', },
    },
    doc => 'Query CG Hub for metadata.',
};

sub help_brief {
}

sub help_deatil {
    return <<HELP;
    https://cghub.ucsc.edu/docs/user/query.html
HELP
}

sub execute {
    my $self = shift;

    my $url = 'https://cghub.ucsc.edu/cghub/metadata/analysisDetail?'.$self->query;
    my $request = HTTP::Request->new(GET => $url);
    my $ua = LWP::UserAgent->new;
    my $response = $ua->request($request);

    die $self->error_message('Failed to execute query! %s', $self->query) if not $response->is_success;

    my $content = $response->content;
    $self->metadata_xml($content);
    Genome::Sys->write_file($self->xml_file, $content) if $self->xml_file;

    return 1;
}

1;

