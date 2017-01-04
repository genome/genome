package Genome::Config::AnalysisProject::Command::TakeOwnership;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::TakeOwnership {
    is => 'Command::V2',
    has_input => [
       analysis_projects  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to take ownership of',
            is_many             => 1,
            shell_args_position => 1,
        }
    ],
};

sub help_brief {
    return 'take ownership of an analysis project';
}

sub help_synopsis {
    return "genome config analysis-project take-ownership <analysis-projects>";
}

sub help_detail {
    return <<"EOS"
Given some analysis projects, this will update the created_by (surrogate for owner) property of the provided analysis projects.
EOS
}

sub execute {
    my $self = shift;
    my $new_user = $ENV{'USER'};
    $self->status_message("Updating created_by user name to: %s", $new_user);
    for my $ap ($self->analysis_projects) {
        my $orig_user = $ap->created_by();
        $ap->created_by($new_user);
        $self->status_message("Updated created_by from %s to %s for project %s", $orig_user, $new_user,$ap->__display_name__);
    }

    return 1;
}

1;
