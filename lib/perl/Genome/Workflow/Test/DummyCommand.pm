package Genome::Workflow::Test::DummyCommand;
class Genome::Workflow::Test::DummyCommand {
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
            default_value => 'apipe',
        },
    ],
};
