package Genome::Ptero::Utils;

use strict;
use warnings FATAL => 'all';
use Params::Validate qw(validate_pos :types);
use Ptero::Proxy::Workflow;
use Try::Tiny qw(try);

use Exporter 'import';

our @EXPORT_OK = qw(
    ptero_proxy
    ptero_workflow_url
);

sub ptero_workflow_url {
    my ($workflow_name) = validate_pos(@_, {type => SCALAR});

    return sprintf("%s?name=%s", Genome::Config::get('ptero_workflow_submit_url'),
        $workflow_name);
}

sub ptero_proxy {
    my ($workflow_name) = validate_pos(@_, {type => SCALAR});

    my $proxy = try {
        return Ptero::Proxy::Workflow->new(
            ptero_workflow_url($workflow_name));
    };
    return $proxy;
}


1;
