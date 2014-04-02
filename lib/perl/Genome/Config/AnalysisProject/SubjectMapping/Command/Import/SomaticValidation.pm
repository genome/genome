package Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation {
    is => 'Command::V2',
    has_input => [
        analysis_project =>  {
            is => 'Genome::Config::AnalysisProject',
            shell_args_position => 1,
            doc => 'The AnalysisProject to associate these experimental pairings with',
        },
        file_path => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'optional path to a newline delimited, tab separated list of sample ids'
        }
    ],
};

sub help_brief {
    return 'Import sample pairings in bulk for an AnalysisProject';
}

sub help_synopsis {
    return "genome config analysis-project subject-mapping import somatic-validation <analysis_project_id> <TSV file path>";
}

sub help_detail {
    return <<"EOS"
This command allows you to import subject mapping information for an AnalysisProject in bulk via a TSV file.

It expects the file to be in a 5 column, tab separated format with the following columns:

tumor_subject_id normal_subject_id snv_variant_list_id indel_variant_list_id sv_variant_list_id

All columns but tumor_subject_id are optional and can be left blank, though the tab separators must be present.
This is useful for setting up single-sample validation models.
A header is optional and should be preceeded with a '#' if present.

EOS
}

my @subjects = ('tumor_subject', 'normal_subject');
my @inputs = ('snv_variant_list_id', 'indel_variant_list_id', 'sv_variant_list_id');

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input                      => $self->file_path,
        separator                  => "\t",
        ignore_lines_starting_with => '#',
        headers                    => [@subjects, @inputs],
    );

    my $count = 0;
    while (my $line = $reader->next()) {
        my $mapping = Genome::Config::AnalysisProject::SubjectMapping->create(analysis_project => $self->analysis_project);

        for(@subjects) {
            my $subject_id = $line->{$_};
            next unless $subject_id;
            $self->_create_subject($mapping, $_, $subject_id);
        }

        for(@inputs) {
            my $value = $line->{$_};
            next unless $value;
            $self->_create_input($mapping, $_, $value);
        }

        $count++;
    }

    $self->status_message("Creating $count new subject mappings.");

    return $count;
}

sub _create_subject {
    my $self = shift;
    my $mapping = shift;
    my $label = shift;
    my $subject_id = shift;

    my $subject = Genome::Subject->get($subject_id);
    die($self->error_message("Unable to find a subject with ID: %s", $subject_id)) unless $subject_id;
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

    my $result = Genome::Model::Tools::DetectVariants2::Result::Base->get($value);
    die($self->error_message("Unable to find a variant list with ID: %s", $value)) unless $result;

    Genome::Config::AnalysisProject::SubjectMapping::Input->create(
        subject_mapping => $mapping,
        value => $value,
        key => $key,
    );
}

1;
