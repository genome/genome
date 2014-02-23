package Genome::Model::Tools::HmpShotgun::CoreMetabolome;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::CoreMetabolome {
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
    'Convert Maq map files into the TCGA format.';
}

sub help_detail {
    return <<EOS
    Convert Maq map files into the TCGA format.
EOS
}


sub execute {
    my $self = shift;
    $self->dump_status_messages(1);
    $self->dump_error_messages(1);
    $self->dump_warning_messages(1);

    my $now = UR::Context->current->now;
    $self->debug_message(">>>Starting CoreMetabalome execute() at $now"); 
    $self->final_file("core_metabolome_final_file_path");
    $self->debug_message("<<<Ending CoreMetabalome execute() at ".UR::Context->current->now); 
    return 1;
}
1;
