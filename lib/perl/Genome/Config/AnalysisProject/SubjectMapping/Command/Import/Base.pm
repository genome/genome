package Genome::Config::AnalysisProject::SubjectMapping::Command::Import::Base;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Import::Base {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        file_path => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'path to a newline-delimited, tab-separated list of subjects (See description section of --help for details.)'
        }
    ],
    is_abstract => 1,
};

sub valid_statuses {
    return ("Pending", "Hold", "In Progress");
}

sub _create_subject {
    my $self = shift;
    my $mapping = shift;
    my $label = shift;
    my $subject_identifier = shift;

    my $subject = Genome::Subject->get($subject_identifier) || Genome::Subject->get(name => $subject_identifier);
    die($self->error_message("Unable to find a subject from identifier: %s", $subject_identifier)) unless $subject;
    Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
        subject_mapping => $mapping,
        subject_id => $subject,
        label => $_,
    );
}

sub _create_input {
    my $self = shift;
    my $mapping = shift;
    my $key = shift;
    my $value = shift;

    Genome::Config::AnalysisProject::SubjectMapping::Input->create(
        subject_mapping => $mapping,
        value => $value,
        key => $key,
    );
}

sub _link_to_tag {
    my $self = shift;
    my $mapping = shift;
    my $tag_string = shift;

    my $tag = Genome::Config::Tag->get($tag_string) || Genome::Config::Tag->get(name => $tag_string);
    die($self->error_message("Unable to find a tag from identifier: %s", $tag_string)) unless $tag;

    $tag->add_subject_mapping($mapping);
}

1;
