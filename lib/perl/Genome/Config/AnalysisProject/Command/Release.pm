package Genome::Config::AnalysisProject::Command::Release;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Release {
    is => 'Command::V2',
    has_input => [
       analysis_projects  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to release for processing',
            is_many             => 1,
            shell_args_position => 1,
        }
    ],
};

sub help_brief {
    return 'resume/start automated model creation for this analysis project';
}

sub help_synopsis {
    return "genome config analysis-project release <analysis-projects>";
}

sub help_detail {
    return <<"EOS"
Given some analysis projects, this will update the status so they can be
automatically processed by CQID.
EOS
}

sub execute {
    my $self = shift;

    for my $ap ($self->analysis_projects) {
        $ap->status('In Progress');
    }

    return 1;
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    for my $status (map{$_->status} $self->analysis_projects){
        unless(grep{$_ eq $status} ("Pending", "Hold", "In Progress")){
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['analysis_projects'],
                desc => "Can't release analysis project with status: $status"
            );
        }
    }
    return @errors;
}

1;
