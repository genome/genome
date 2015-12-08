package Genome::Config::AnalysisProject::SubjectMapping::Command::Export::SomaticValidation;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Export::SomaticValidation {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        file_path => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'path to write the output file',
            default_value => '-',
        },
        include_header => {
            is => 'Boolean',
            doc => 'include a header line in the file',
            default_value => 0,
            is_optional => 1,
        },
    ],
};

sub help_brief {
    return 'Export sample pairings in bulk from an AnalysisProject';
}

sub help_synopsis {
    return "genome config analysis-project subject-mapping export somatic-validation <analysis_project_id> <TSV file path>";
}

sub help_detail {
    return <<"EOS"
This command allows you to export subject mapping information from an AnalysisProject in bulk via a TSV file.

It produces a file in a 5+ column, tab separated format with the following columns:

tumor_subject normal_subject snv_variant_list_id indel_variant_list_id sv_variant_list_id [tag...]

Only a tumor_subject is required, but the tab separators will be present as placeholders for the normal_subject and variant lists.  If any number of tags are present, they will be listed as extra columns beginning with the sixth.
EOS
}

my @column_order = (qw(tumor_sample normal_sample snv_variant_list_id indel_variant_list_id sv_variant_list_id));

sub valid_statuses {
    #Any status is fine for export.
    return @{Genome::Config::AnalysisProject->__meta__->property(property_name => 'status')->valid_values};
}

sub execute {
    my $self = shift;

    my @subject_mappings = Genome::Config::AnalysisProject::SubjectMapping->get(analysis_project => $self->analysis_project);

    my $fh = Genome::Sys->open_file_for_writing($self->file_path);
    $self->write_line($fh, '#'.$column_order[0], @column_order[1..$#column_order]) if $self->include_header;

    for my $sm (@subject_mappings) {
        my %data;

        my @subjects = $sm->subject_bridges;
        for my $s (@subjects) {
            $data{$s->label} = $s->subject->id;
        }

        my @inputs = $sm->inputs;
        for my $i (@inputs) {
            $data{$i->key} = $i->value;
        }

        my @tags = $sm->tags;
        if(@tags) {
            $data{tags} = join("\t", map($_->name, @tags));
        }

        $self->write_line($fh, @data{@column_order}, ($data{tags}? $data{tags}: ()));
    }

    return 1;
}

sub write_line {
    my $self = shift;
    my $fh = shift;
    my @values = @_;

    $fh->say(join("\t", @values));
}

1;
