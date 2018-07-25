package Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation {
    is => 'Genome::Config::AnalysisProject::SubjectMapping::Command::Import::Base',
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

It expects the file to be in a 5+ column, tab separated format with the following columns:

tumor_subject normal_subject snv_variant_list_id indel_variant_list_id sv_variant_list_id [tag...]

Only a tumor_subject is required, but the tab separators must be present as placeholders for the normal_subject and variant lists.  If any number of tags are desired, they can be listed as extra columns beginning with the sixth.
A header is optional and should be preceded with a '#' if present.
Both tumor and normal subject can be specified by either ID or Name.
Tags may also be specified by either ID or Name.
EOS
}

my @subjects = ('tumor_sample', 'normal_sample');
my @inputs = ('snv_variant_list_id', 'indel_variant_list_id', 'sv_variant_list_id');

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input                      => $self->file_path,
        separator                  => "\t",
        ignore_lines_starting_with => '#',
        headers                    => [@subjects, @inputs],
        allow_extra_columns        => 1,
    );

    my $count = 0;
    while (my $line = $reader->next()) {
        my $mapping = Genome::Config::AnalysisProject::SubjectMapping->create(analysis_project => $self->analysis_project);

        for(@subjects) {
            my $subject_identifier = $line->{$_};
            next unless $subject_identifier;
            $self->_create_subject($mapping, $_, $subject_identifier);
        }

        for(@inputs) {
            my $value = $line->{$_};
            next unless $value;
            $self->_create_input($mapping, $_, $value);
        }

        my $tags = $reader->current_extra_columns;
        for(@$tags) {
            $self->_link_to_tag($mapping, $_);
        }

        $count++;
    }

    $self->status_message("Creating $count new subject mappings.");

    return $count;
}


sub _create_input {
    my $self = shift;
    my $mapping = shift;
    my $key = shift;
    my $value = shift;

    my $result = Genome::Model::Tools::DetectVariants2::Result::Base->get($value);
    die($self->error_message("Unable to find a variant list with ID: %s", $value)) unless $result;

    $self->SUPER::_create_input($mapping, $key, $value);
}


1;
