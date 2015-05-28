package Genome::WorkflowBuilder::Test::DummyCommandCalculatedDefault;
class Genome::WorkflowBuilder::Test::DummyCommandCalculatedDefault {
    is => [qw(Command Genome::Configurable)],

    has_input => ['input'],

    has_output => [
        many_output => {
            is_many => 1,
        },
        single_output => { },
    ],

    has_param => [
        lsf_resource => {
            config => 'dummy_lsf_resource',
        },
        lsf_queue => {
            config => 'dummy_lsf_queue',
        },
    ],
};
