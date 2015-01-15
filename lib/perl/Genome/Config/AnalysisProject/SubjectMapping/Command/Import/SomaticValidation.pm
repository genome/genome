package Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Command::Import::SomaticValidation {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        file_path => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'path to a newline-delimited, tab-separated list of samples, variant lists, and tags (See description section of --help for details.)'
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

It expects the file to be in a 5+ column, tab separated format with the following columns:

tumor_subject normal_subject snv_variant_list_id indel_variant_list_id sv_variant_list_id [tag...]

Only a tumor_subject is required, but the tab separators must be present as placeholders for the normal_subject and variant lists.  If any number of tags are desired, they can be listed as extra columns beginning with the sixth.
A header is optional and should be preceded with a '#' if present.
Both tumor and normal subject can be specified by either ID or Name.
Tags may also be specified by either ID or Name.
EOS
}

sub valid_statuses {
    return ("Pending", "Hold", "In Progress");
}

my @subject_labels = ('tumor_sample', 'normal_sample');
my @inputs = ('snv_variant_list_id', 'indel_variant_list_id', 'sv_variant_list_id');

sub execute {
    my $self = shift;

    my $reader = Genome::Utility::IO::SeparatedValueReader->create(
        input                      => $self->file_path,
        separator                  => "\t",
        ignore_lines_starting_with => '#',
        headers                    => [@subject_labels, @inputs],
        allow_extra_columns        => 1,
    );

    my $count = 0;
    while (my $line = $reader->next()) {
        my $mapping = Genome::Config::AnalysisProject::SubjectMapping->create(analysis_project => $self->analysis_project);

        $self->_add_subject_bridges($mapping, $line);

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

sub _add_subject_bridges {
    my ($self, $mapping, $line) =@_;

    if ( $line->{normal_sample} eq $line->{tumor_sample} ) {
        die $self->error_message('Same sample given for normal/tumor: %s', $line->{normal_sample});
    }

    my %subject_common_names_seen;
    for my $subject_label (@subject_labels) {
        my $subject_identifier = $line->{$subject_label};
        next unless $subject_identifier;
        my $subject_common_name = $self->_create_subject($mapping, $subject_label, $subject_identifier);
        $subject_common_names_seen{$subject_common_name.'_sample'}++;
    }

    return 1;
}

sub _create_subject {
    my $self = shift;
    my $mapping = shift;
    my $label = shift;
    my $subject_identifier = shift;

    my $subject = Genome::Subject->get($subject_identifier) || Genome::Subject->get(name => $subject_identifier);
    die($self->error_message("Unable to find a subject from identifier: %s", $subject_identifier)) unless $subject;

    my $subject_common_name = $self->determine_if_subject_is_norml_or_tumor($subject);
    if ( $subject_common_name eq 'normal' ) { # only checking if normal sample is labeled tumor
        my $expected_label = $subject_common_name.'_sample';
        if ( $label ne $expected_label ) {
            die $self->error_message('Incompatible label (%s) for subject %s!', $label, $subject->__display_name__);
        }
    }

    my $subject_bridge = Genome::Config::AnalysisProject::SubjectMapping::Subject->create(
        subject_mapping => $mapping,
        subject_id => $subject,
        label => $_,
    );
    die $self->error_message('Failed to create %s subject bridge for subject (%s)!', $label, $subject->__display_name__) if not $subject_bridge;

    return $subject_common_name;
}

sub determine_if_subject_is_norml_or_tumor {
    my ($self, $subject) = @_;

    die 'No subject given!' if not $subject;

    my $common_name = $subject->common_name;
    return 'unknown' if not defined $common_name;
    return 'normal' if $common_name =~ /normal/i;
    for my $possible_tumor_common_name (qw/ tumor xenograft sarcoma relapse /) {
        return 'tumor' if $common_name =~ /$possible_tumor_common_name/i;
    }

    return 'unknown';
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

sub _link_to_tag {
    my $self = shift;
    my $mapping = shift;
    my $tag_string = shift;

    my $tag = Genome::Config::Tag->get($tag_string) || Genome::Config::Tag->get(name => $tag_string);
    die($self->error_message("Unable to find a tag from identifier: %s", $tag_string)) unless $tag;

    $tag->add_subject_mapping($mapping);
}

1;
