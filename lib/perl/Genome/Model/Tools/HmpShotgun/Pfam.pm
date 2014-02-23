package Genome::Model::Tools::HmpShotgun::Pfam;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::Pfam {
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
    'Pfam.';
}

sub help_detail {
    return <<EOS
    Pfam.
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

    my $now = UR::Context->current->now;
    $self->debug_message(">>>Starting Pfam execute() at $now"); 
    $self->final_file("pfam_final_file_path");
    $self->debug_message("<<<Ending Pfam execute() at ".UR::Context->current->now); 
    return 1;

}
1;
