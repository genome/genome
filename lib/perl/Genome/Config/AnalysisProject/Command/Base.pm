package Genome::Config::AnalysisProject::Command::Base;

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Config::AnalysisProject::Command::Base {
    is => 'Command::V2',
    is_abstract => 1,
    #has => [ analysis_project... ], #added in preprocessor
    has_transient => [
        valid_statuses => {
            is => 'ARRAY',
            doc => 'List of analysis project statuses valid for running this command',
            default_value => [],
        },
    ],
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
};

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    unless (exists $desc->{has}{valid_statuses} and $desc->{has}{valid_statuses}{default_value}) {
        die 'Bad class definition for ' . $class . '.  Commands inheriting from Genome::Config::AnalysisProject::Command::Base must define a "valid_statuses" attribute.';
    }

    my $statuses = $desc->{has}{valid_statuses}{default_value};

    $desc->{has}{analysis_project} = {
        is => 'Genome::Config::AnalysisProject',
        shell_args_position => 1,
        doc => 'the analysis project on which to operate--must be in one of the following statuses: ' . join(', ', @$statuses),
    };

    return $desc;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);

    my $status = $self->analysis_project->status;
    unless(grep{$_ eq $status} @{$self->valid_statuses}){
        my $name = $self->command_name;
        $name = (split(' ', $name))[-1];
        $name =~ s/-/ /g;

        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['analysis_project'],
            desc => sprintf(q(Can't %s using analysis project with status: %s), $name, $status),
        );
    }

    return @errors;
}

1;
