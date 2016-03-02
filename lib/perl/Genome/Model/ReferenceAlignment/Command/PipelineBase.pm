package Genome::Model::ReferenceAlignment::Command::PipelineBase;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Command::PipelineBase {
    is => 'Command::V2',
    is_abstract => 1,
    has_input_output => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            doc => 'build for which to run alignments',
        },
    ],
    doc => 'base class for reference alignment pipeline commands',
};

sub sub_command_category { 'pipeline steps' }

1;
