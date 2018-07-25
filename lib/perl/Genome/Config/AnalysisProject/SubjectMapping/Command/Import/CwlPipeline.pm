package Genome::Config::AnalysisProject::SubjectMapping::Command::Import::CwlPipeline;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Import::CwlPipeline {
    is => 'Genome::Config::AnalysisProject::SubjectMapping::Command::Import::Base',
};

sub help_brief {
    return 'Import subject mappings in bulk for an AnalysisProject';
}

sub help_synopsis {
    return "genome config analysis-project subject-mapping import cwl-pipeline <analysis_project_id> <TSV file path>";
}

sub help_detail {
    return <<"EOS"
This command allows you to import subject mapping information for an AnalysisProject in bulk via a TSV file.

Each line of the file should provide information about a single subject-mapping.

For each subject in the mapping, include two columns: the first with the label and the second with the identifier of the subject.
After the subject(s), optionally include inputs as additional two-column pairs: the first with the "key" prefixed by "i:" and the second with the value.
Additionally, optionally include any tags as additional columns with the tag prefixed by "t:".

Example line:

tumor_sample\t1234567890\tnormal_sample\t9876543210\ti:prior_variants_result\tabcdef1234567890fedcba9876543210\tt:discovery

Only a single subject is required.  Each line can have a variable number of subjects and inputs.
A header is optional and should be preceded with a '#' if present.
Subject can be specified by either ID or Name (but using a Name may fail if multiple subjects share it, such as when the Individual and Sample have the same Name.)
Tags may also be specified by either ID or Name.
EOS
}

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input                      => $self->file_path,
        separator                  => "\t",
        ignore_lines_starting_with => '#',
        headers                    => [],
        allow_extra_columns        => 1,
    );

    my $count = 0;
    while ($reader->next()) {
        my $mapping = Genome::Config::AnalysisProject::SubjectMapping->create(analysis_project => $self->analysis_project);

        #since we passed no headers, all of the data is here
        my $data = $reader->current_extra_columns;

        while (@$data) {
            my $next = shift @$data;

            if ($next =~ s/^([it])://) {
                my $type = $1;
                if ($type eq 'i') {
                    if (@$data) {
                        my $value = shift @$data;
                        $self->_create_input($mapping, $next, $value);
                    } else {
                        $self->fatal_message('Missing value for input specification %s on line %s', $next, $reader->line_number);
                    }
                } elsif ($type eq 't') {
                    $self->_link_to_tag($mapping, $next);
                } else {
                    die 'program error: matched something unexpected in regex!';
                }
            } else {
                if (@$data) {
                    my $identifier = shift @$data;
                    $self->_create_subject($mapping, $next, $identifier);
                } else {
                    $self->fatal_message('Missing identifier for subject specification %s on line %s', $next, $reader->line_number);
                }
            }
        }

        $count++;
    }

    $self->status_message('Creating %s new subject mappings.', $count);

    return $count;
}

1;
