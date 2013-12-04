package Genome::WorkflowBuilder::Test::DummyCommand;
class Genome::WorkflowBuilder::Test::DummyCommand {
    is => 'Command',

    has_input => ['input'],

    has_output => [
        many_output => {
            is_many => 1,
        },
        single_output => { },
    ],

    has_param => [
        lsf_resource => {
            default_value => "-M 25000000 -R 'select[mem>25000] rusage[mem=25000]'",
        },
        lsf_queue => {
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
    ],
};
