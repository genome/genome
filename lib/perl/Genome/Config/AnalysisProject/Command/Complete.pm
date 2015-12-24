package Genome::Config::AnalysisProject::Command::Complete;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Complete {
    is => 'Command::V2',
    has_input => [
       analysis_projects  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis projects to complete',
            is_many             => 1,
            shell_args_position => 1,
        }
    ],
};

sub help_brief {
    return 'complete this analysis project';
}

sub help_synopsis {
    return "genome config analysis-project complete <analysis-projects>";
}

sub help_detail {
    return <<"EOS"
Given some analysis projects, this will update the status so CQID will no longer produce new models.
EOS
}

sub execute {
    my $self = shift;

    for my $ap ($self->analysis_projects) {
        $ap->status('Completed');
        $self->status_message("Completed: %s", $ap->__display_name__);
    }

    return 1;
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    for my $status (map{$_->status} $self->analysis_projects){
        unless(grep{$_ eq $status} ("In Progress")){
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => ['analysis_projects'],
                desc => "Can't complete analysis project with status: $status"
            );
        }
    }
    return @errors;
}

1;
