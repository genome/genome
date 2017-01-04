package Genome::Config::AnalysisProject::Command::Uncomplete;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Uncomplete {
    is => 'Command::V2',
    has_input => [
       analysis_projects  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to complete',
            is_many             => 1,
            shell_args_position => 1,
            require_user_verify => 1,
        }
    ],
};

sub help_brief {
    return 'uncomplete this analysis project';
}

sub help_synopsis {
    return "genome config analysis-project uncomplete <analysis-projects>";
}

sub help_detail {
    return <<"EOS"
Given some analysis projects you created, this will revert them to "In Progress" from a completed state.
EOS
}

sub execute {
    my $self = shift;

    for my $ap ($self->analysis_projects) {
        $ap->status('In Progress');
        $self->status_message("Reverted AnP to 'In Progress': %s", $ap->__display_name__);
    }

    return 1;
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    for my $anp ($self->analysis_projects) {
        unless($anp->status eq "Completed"){
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['analysis_projects'],
                desc => 'Cannot process incomplete analysis project ' . $anp->id,
            );
        }
        unless($anp->created_by eq Genome::Sys->username){
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['analysis_projects'],
                desc => 'Cannot process analysis project created by someone else: ' . $anp->id,
            );
        }
    }
    return @errors;
}

1;
