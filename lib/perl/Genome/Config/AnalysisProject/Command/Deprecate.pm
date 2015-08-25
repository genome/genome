package Genome::Config::AnalysisProject::Command::Deprecate;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Deprecate {
    is => 'Genome::Config::AnalysisProject::Command::Base',
};

sub help_brief {
    return 'deprecate the analysis project';
}

sub help_synopsis {
    return "genome config analysis-project deprecate <analysis-project>";
}

sub help_detail {
    return <<"EOS"
Deprecate the input analysis project.  Deprecate includes setting all config profile items to disabled and updating the project status to Deprecated.
EOS
}

sub valid_statuses {
    return ("Pending", "Hold", "In Progress", "Template");
}

sub execute {
    my $self = shift;

    my $anp = $self->analysis_project;
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

    return 1;
}


1;
