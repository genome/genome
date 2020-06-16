package Genome::Model::CwlPipeline::Command::ListDockerImages;

use strict;
use warnings;

use Genome;

use Cwd qw();
use File::Basename qw();
use File::Spec;
use JSON qw();

class Genome::Model::CwlPipeline::Command::ListDockerImages {
    is => 'Command::V2',
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile::CwlPipeline',
            doc => 'The processing profile for which to collect information',
            shell_args_position => 1,
        },
    ],
    has_optional => [
        output_file => {
            is => 'File',
            doc => 'where to write the output (or STDOUT if not specified',
        },
    ],
    has_transient_optional => {
        _fh => {
            doc => 'the filehandle for output',
        },
    },
    doc => 'list the docker images used by a CWL Pipeline',
};

sub help_detail {
    return <<EOHELP
Generates a list of all the docker images specified in a given workflow.
EOHELP
    ;
}

sub execute {
    my $self = shift;

    my $out = $self->output_file // '-';
    $self->_fh(Genome::Sys->open_file_for_writing($out));

    my $pp = $self->processing_profile;
    my $main_workflow_file = $pp->main_workflow_file;

    unless ($main_workflow_file =~ /\.cwl$/) {
        $self->fatal_message('Only CWL pipelines are currently supported by this tool.');
    }

    $self->process_file($main_workflow_file);

    $self->_fh->close if $self->output_file;

    return 1;
}

sub process_file {
    my $self = shift;
    my $file = shift;

    my $canonical_cwl = `/usr/local/bin/cwltool --pack $file`;
    my $data = JSON::decode_json($canonical_cwl);

    my $steps = $data->{'$graph'};
    unless ($steps and ref($steps) eq 'ARRAY' and @$steps) {
        $self->fatal_message(q{Couldn't get canonical representation of CWL for file: %s}, $file);
    }

    for my $step (@$steps) {
        $self->process_step($step);
    }

    return 1;
}

sub process_step {
    my $self = shift;
    my $data = shift;

    my $id = $data->{id};
    $self->output_id($id);

    if (my $reqs = $data->{requirements}) {
        $self->process_requirements($reqs);
    }
    if (my $hints = $data->{hints}) {
        $self->process_requirements($hints);
    }

    return 1;
}

sub process_requirements {
    my $self = shift;
    my $requirements = shift;

    return unless $requirements and @$requirements;

    for my $req (@$requirements) {
        next unless $req->{class} eq 'DockerRequirement';

        if (my $image = $req->{dockerImageId}) {
            $self->output_image($image);
        }

        if (my $image = $req->{dockerPull}) {
            $self->output_image($image);
        }
    }

    return 1;
}

sub output_id {
    my $self = shift;
    my $id = shift;

    $self->_fh->say($id);
}

sub output_image {
    my $self = shift;
    my $image = shift;

    $self->_fh->say('  ', $image);
}

1;
