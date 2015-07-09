package Genome::Config::AnalysisProject::Command::Deprecate;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Deprecate {
    is => 'Command::V2',
    has_input => [
       analysis_projects  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to deprecate',
            is_many             => 1,
            shell_args_position => 1,
        }
    ],
};

sub help_brief {
    return 'deprecate the analysis project';
}

sub help_synopsis {
    return "genome config analysis-project deprecate <analysis-projects>";
}

sub help_detail {
    return <<"EOS"
Given some analysis projects, this will deprecate the analysis project.  This includes setting all config profile items to disabled and updating the project status to Deprecated.
EOS
}

sub execute {
    my $self = shift;

    for my $anp ($self->analysis_projects) {
        my @active_config_profile_items = grep {$_->status eq 'active'} $anp->config_items;
        for my $config_profile_item (@active_config_profile_items) {
            my $cmd = Genome::Config::AnalysisProject::Command::DisableConfigFile->create(
                profile_item => $config_profile_item,
            );
            unless ( $cmd->execute() ) {
                $self->error_message('Failed to disable config profile item: '. $config_profile_item->id);
                die($self->error_message);
            }
        }
        $anp->status('Deprecated');
    }

    return 1;
}


1;
