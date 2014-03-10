package Genome::Model::Tools::HmpShotgun::UnalignedReadAnalysis;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::UnalignedReadAnalysis {
    is  => ['Command'],
    has => [
        model_id => {
            is  => 'String',
            is_input => '1',
            doc => 'The model id to process.',
        },
        working_directory => {
            is => 'String',
            is_input => '1',
            doc => 'The working directory where results will be deposited.',
        },
         final_file => {
            is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The model id to process.',
        },
        delete_intermediates => {
            is => 'Integer',
            is_input =>1,
            is_optional =>1,
            default=>0,
        },
    ],
};


sub help_brief {
    'Handle the unaligned reads in a novel fashion.';
}

sub help_detail {
    return <<EOS
    Handle the unaligned reads in a novel fashion.
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

    my $now = UR::Context->current->now;
    $self->debug_message(">>>Starting Unaligned execute() at $now"); 
    $self->final_file("unaligned_final_file_path");
    $self->debug_message("<<<Ending Unaligned execute() at ".UR::Context->current->now); 
    return 1;

}
1;
