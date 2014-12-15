package Genome::Config::AnalysisProject::Command::Base;

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Config::AnalysisProject::Command::Base {
    is => 'Command::V2',
    is_abstract => 1,
    #has => [ analysis_project... ], #added in preprocessor
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
};

sub valid_statuses {
    my $class = shift;
    $class = ref $class || $class;

    die sprintf('Class %s must define a "valid_statuses" method indicating those statuses an analysis project can be in to be used by the command.', $class);
}

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;

    my @statuses = $desc->{class_name}->valid_statuses;

    $desc->{has}{analysis_project} = {
        is => 'Genome::Config::AnalysisProject',
        shell_args_position => 1,
        doc => 'the analysis project on which to operate--must be in one of the following statuses: ' . join(', ', @statuses),
    };

    return $desc;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);

    my $status = $self->analysis_project->status;
    unless(grep{$_ eq $status} $self->valid_statuses){
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
