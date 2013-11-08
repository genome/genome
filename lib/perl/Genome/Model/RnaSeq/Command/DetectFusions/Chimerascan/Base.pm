package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::Base {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Base',
};

sub command_class_prefix {
    die "Abstract";
}

sub execute {
    my $self = shift;

    my $dag = $self->build_workflow();
    print Data::Dumper::Dumper($dag->get_xml());
    my $outputs = $dag->execute(%{$self->workflow_inputs});

    $self->result($outputs->{'result'});
    return 1;
}

sub workflow_inputs {
    my $self = shift;

    my ($bowtie_version, $reuse_bam, $detector_params) = $self->_resolve_options;

    my %inputs = (
        detector_version => $self->detector_version,
        detector_params => $detector_params,
        bowtie_version => $bowtie_version,
        reuse_bam => $reuse_bam,
        build => $self->build,
    );
    return \%inputs;
}

sub build_workflow {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'ChimerascanWorkflow',
        log_dir => $self->build->log_directory,
    );

    my $index_command = $self->attach_index_command_to($dag);
    my $detector_command = $self->attach_detector_command_to($dag);

    $self->_add_links($dag, $index_command, $detector_command);
    $self->_add_outputs($dag, $detector_command);

    return $dag;
}

sub attach_detector_command_to {
    my $self = shift;
    my $dag = shift;

    my $detector_command = Genome::WorkflowBuilder::Command->create(
        name => 'DetectorCommand',
        command => $self->detector_command_name,
    );
    $dag->add_operation($detector_command);
    $self->_add_common_inputs($dag, $detector_command);
    $dag->connect_input(
        input_property => 'reuse_bam',
        destination => $detector_command,
        destination_property => 'reuse_bam',
    );
    return $detector_command;
}

sub attach_index_command_to {
    my $self = shift;
    my $dag = shift;

    my $index_command = Genome::WorkflowBuilder::Command->create(
        name => 'IndexCommand',
        command => $self->index_command_name,
    );
    $dag->add_operation($index_command);
    $self->_add_common_inputs($dag, $index_command);
    $dag->connect_input(
        input_property => 'build',
        destination => $index_command,
        destination_property => 'build',
    );
    return $index_command;
}

sub _add_common_inputs {
    my $self = shift;
    my $dag = shift;
    my $command = shift;

    my @common_inputs = qw(
        detector_version
        detector_params
        bowtie_version
    );

    for my $prop_name (@common_inputs) {
        $dag->connect_input(
            input_property => $prop_name,
            destination => $command,
            destination_property => $prop_name,
        );
    }
}

sub _add_links {
    my $self = shift;
    my $dag = shift;
    my $source_op = shift;
    my $destination_op = shift;

    $dag->create_link(
        source => $source_op, source_property => 'build',
        destination => $destination_op, destination_property => 'build',
    );
}

sub _add_outputs {
    my $self = shift;
    my $dag = shift;
    my $command = shift;

    $dag->connect_output(
        output_property => 'result',
        source => $command,
        source_property => 'result',
    );
}

sub detector_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::Detector";
}

sub index_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::Index";
}

# these are the options which you must specify to us in the
# fusion-detection-strategy part of the processing-profile
our %OUR_OPTIONS_VALIDATORS = (
        '--bowtie-version' => '_validate_bowtie_version',
        '--reuse-bam' => '_validate_reuse_bam',
);

sub _validate_bowtie_version {
    my ($self, $val, $params) = @_;

    my $bowtie_version = $val ||
            die("You must supply a bowtie version in the detector parameters " .
                "in the form of \"--bowtie-version=<version>\" Got detector " .
                "parameters: [\"$params\"]");
    my ($major_version) = split(/\./, $bowtie_version);
    if ($major_version ne 0) {
        die("Currently chimerascan only supports bowtie major version 0, " .
            "not $major_version");
    }
}

sub _validate_reuse_bam {
    my ($self, $val, $params) = @_;

    unless ($val eq 1 or $val eq 0) {
        die "You must specify either 1 (true) or 0 (false) for parameter " .
                \"--reuse-bam\", you specified \"$val\"";
    }
}


# return our options (hash) and the options for chimerascan (string)
sub _resolve_options {
    my $self = shift;

    my $params = $self->detector_params;
    my %our_opts;
    # go through and remove our options from the params
    for my $name (keys %OUR_OPTIONS_VALIDATORS) {
        # \Q$foo\E ensures that regex symbols are 'quoted'
        if($params and $params =~ m/(\s*\Q$name\E[=\s]([^\s]*)\s*)/) {
            my $str = $1;
            my $val = $2;
            $params =~ s/\Q$str\E/ /;

            my $validation_method_name = $OUR_OPTIONS_VALIDATORS{$name};
            $self->$validation_method_name($val, $params); # dies if invalid

            $our_opts{$name} = $val;
        } else {
            my $t = q(Could not find parameter named '%s' in param string '%s');
            die(sprintf($t, $name, $params || ''));
        }
    }
    return ($our_opts{'--bowtie-version'}, $our_opts{'--reuse-bam'}, $params);
}

1;
