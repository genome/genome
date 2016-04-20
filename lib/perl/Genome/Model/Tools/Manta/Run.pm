package Genome::Model::Tools::Manta::Run;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Manta::Run {
    is => 'Genome::Model::Tools::Manta::Base',
    has_input => [
        working_directory => {
            is => 'Text',
            doc => 'The same working directory used during the Config step.  The directory should contain the  script used to execute Manta',
        },
    ],
    has_optional_param => [
        mode => {
            is => 'Text',
            doc => 'select run mode',
            tool_arg_name => 'mode',
            # default values are not optimal, but we only want to run in local mode since sge is not a valid option in the GMS
            default_value => 'local',
            is_mutable => 0,
        },
        memory => {
            is => 'Number',
            doc => 'gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer ( Estimate the total memory for this node for local mode).',
            tool_arg_name => 'memGb',
        },
        email => {
            is => 'Text',
            doc => 'send email notification of job completion status to this address (may be provided multiple times for more than one email address)',
            is_many => 1,
            tool_arg_name => 'mailTo',
        },
        dry_run => {
            is => 'Boolean',
            doc => 'dryRun workflow code without actually running command-tasks',
            tool_arg_name => 'dryRun',
        },
        quiet => {
            is => 'Boolean',
            doc => 'Do not write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)',
            tool_arg_name => 'quiet',
        },
        # The version is not required since the python paths are configured separately
        version => {
        },
        # Ignoring manta params for SGE : queue and jobs
    ],
};

# Specifically over-riding the Base class implementation to execute python script in working directory
sub tool_path {
    my $self = shift;
    return File::Spec->join($self->working_directory,$self->_tool_subcommand_name);
}

sub _tool_subcommand_name {
    return 'runWorkflow.py';
}

1;
