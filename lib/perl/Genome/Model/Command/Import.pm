package Genome::Model::Command::Import;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Import {
    is => 'Genome::Model::Command::BaseDeprecated',
    has => [ 
           model => { is => 'Genome::Model', id_by => 'model_id'},
       ],
};

sub help_brief {
    "export data about a model into an external format"
}
sub generate_report_brief {
    die "Implement me in the subclass, owangutang";
}
sub generate_report_detail {
    die "Implement generate_report_detail in the subclass you're writing, or this will 
     continue to fail";
}

sub resolve_reports_directory {
    my $self = shift;
    my $model = $self->model;
    my $last_build = $model->last_complete_build;
    unless ($last_build) {
        $self->error_message("No successful build of " . $model->id);
        return;
    }

    my $reports_dir = $last_build->resolve_reports_directory;
    unless (-d $reports_dir) {
        unless (mkdir $reports_dir) {
            $self->error_message("Directory $reports_dir doesn't exist, can't create");
            return;
        }
        chmod 02770, $reports_dir;
    }

    return $reports_dir;
}

1;
