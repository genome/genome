package Genome::Config::AnalysisProject::Command::TakeOwnership;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::TakeOwnership {
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
    return 'take analysis projects';
}

sub help_synopsis {
    return "genome config analysis-project take <analysis-projects>";
}

sub help_detail {
    return <<"EOS"
Given some analysis projects, this will update the "created_by" to the current user.
EOS
}

sub execute {
    my $self = shift;

    for my $ap ($self->analysis_projects) {
        $ap->created_by(Genome::Sys->username);
    }

    return 1;
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    for my $anp ($self->analysis_projects){
        if ($anp->is_cle) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['analysis_projects'],
                desc => 'Cannot take CLE analysis project ' . $anp->id,
            );
        }

        if (Genome::Sys->user_has_role($anp->created_by, 'production')) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['analysis_projects'],
                desc => 'Cannot take production analysis project ' . $anp->id,
            );
        }
    }
    return @errors;
}

1;
