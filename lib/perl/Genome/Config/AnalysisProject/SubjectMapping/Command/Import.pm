package Genome::Config::AnalysisProject::SubjectMapping::Command::Import;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Import {
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
    return "genome analysis-project subject-mapping import <analysis_project_id> <TSV file path>";
}

sub help_detail {
    return <<"EOS"
This command allows you to import subject mapping information for an AnalysisProject in bulk.

It expects each line of the file to be formatted as follows:
control_subject_id<TAB>experimental_subject_id
EOS
}

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input     => $self->file_path,
        headers   => ['control_subject_id', 'experimental_subject_id'],
        separator => "\t",
    );

    my $count = 0;
    while (my $line = $reader->next()) {
        Genome::Config::AnalysisProject::SubjectMapping->create(
            %$line,
            analysis_project => $self->analysis_project,
        );
        $count++;
    }

    $self->status_message("Creating $count new subject mappings.");

    return $count;
}

1;
