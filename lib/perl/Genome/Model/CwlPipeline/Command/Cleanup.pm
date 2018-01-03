package Genome::Model::CwlPipeline::Command::Cleanup;

use strict;
use warnings;

use File::Spec;
use Genome;

class Genome::Model::CwlPipeline::Command::Cleanup {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build::CwlPipeline',
            is_many => 1,
            doc => 'The build(s) to cleanup disk.',
        },
    ],
    doc => 'cleanup disk for older CwlPipeline builds.',
};

sub help_detail {
    return <<EOHELP
Remove the tmp* directories and the content of those directories from the build data directory. Also reallocate the size after removal.
EOHELP
    ;
}

sub execute {
    my $self = shift;

    for my $build ($self->builds) {
        unless ($build->status eq 'Succeeded') {
            $self->fatal_message("Unable to run cleanup on '%s' build with id '%s'. For Failed builds, please abandon instead once troubleshooting is complete.",$build->status,$build->id);
        }
        
        my $data_directory = $build->data_directory;

        my $results_dir = File::Spec->join($data_directory, 'results');

        my $tmp_dir = Genome::Sys->create_temp_directory($build->id);
        unless (Genome::Model::CwlPipeline::Command::Run->cleanup($tmp_dir, $results_dir)) {
            $self->fatal_message("Failed to cleanup build tmp dir '%s' and results dir '%s'", $tmp_dir, $results_dir);
        }

        my $allocation_path = File::Spec->join('model_data',$build->model->id,'build'. $build->id);
        my $allocation = Genome::Disk::Allocation->get(allocation_path => $allocation_path);
        unless ($allocation) {
            $self->fatal_message("Failed to find allocation by allocation path '%s'.", $allocation_path);
        }
        $allocation->reallocate();
    }

    return 1;
}


1;
