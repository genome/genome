package Genome::Model::Build::Msi;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Msi {
    is => 'Genome::Model::Build',
};

sub calculate_estimated_kb_usage {
    my $self = shift;
    
    my $assembly_model = $self->model->from_models;

    my $assembly_build = $assembly_model->last_complete_build;
    unless ($assembly_build) {
        $self->error_message("Underlying model " . $assembly_model->__display_name__ . " has no complete builds!");
        return;
    }
    
    my $assembly_build_directory = $assembly_build->data_directory;
    unless ($assembly_build_directory && -e $assembly_build_directory) {
        my $msg = $self->error_message("Failed to get last complete build directory for the input assembly!");
        Carp::confess($msg);
    }
    
    my $du_out = `du -sk $assembly_build_directory/edit_dir`;
    my $input_build_kb = 512_000;
    ($input_build_kb) = $du_out =~ /(\d+).*/;
    
    return POSIX::floor($input_build_kb*2.5);
}
1;



