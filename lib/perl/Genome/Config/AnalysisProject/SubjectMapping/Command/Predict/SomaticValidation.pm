package Genome::Config::AnalysisProject::SubjectMapping::Command::Predict::SomaticValidation;

use strict;
use warnings;

use Genome;
use List::MoreUtils qw(uniq);

class Genome::Config::AnalysisProject::SubjectMapping::Command::Predict::SomaticValidation {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        file_path => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'output path for a newline-delimited, tab-separated list of samples and required whitespace'
        }
    ],
    has_optional => [
        normal_sample_common_names => {
            is      => 'Text',
            default => Genome::Model::ClinSeq::Command::UpdateAnalysis->__meta__->property('normal_sample_common_names')->default_value,
            doc     => Genome::Model::ClinSeq::Command::UpdateAnalysis->__meta__->property('normal_sample_common_names')->doc,
        },
        tumor_sample_common_names => {
            is      => 'Text',
            default => Genome::Model::ClinSeq::Command::UpdateAnalysis->__meta__->property('tumor_sample_common_names')->default_value,
            doc => Genome::Model::ClinSeq::Command::UpdateAnalysis->__meta__->property('tumor_sample_common_names')->doc,
        },
        sample_identifier => {
            is => 'Text',
            doc => 'The method used to resolve the identifier of the sample.',
            default => 'name',
            valid_values => ['name','id'],
        },
    ],
};

sub help_brief {
    return 'Predict the sample pairings needed to create Somatic Validation models for an AnalysisProject';
}

sub help_synopsis {
    return "genome config analysis-project subject-mapping predict somatic-validation <analysis_project_id> <TSV file path>";
}

sub help_detail {
    return <<"EOS"
This command will predict subject mapping for Somatic Validation model types across all samples in an AnalysisProject.

The output file is a 6 column, tab separated format with the following columns:

tumor_subject normal_subject snv_variant_list_id indel_variant_list_id sv_variant_list_id tag

EOS
}

sub valid_statuses {
    return ("Pending", "Hold", "In Progress");
}

my @headers = ('tumor_sample', 'normal_sample','snv_variant_list_id', 'indel_variant_list_id', 'sv_variant_list_id','tag');

sub execute {
    my $self = shift;

    my $analysis_project = $self->analysis_project;

    my $normal_sample_common_names = $self->normal_sample_common_names;
    my $tumor_sample_common_names  = $self->tumor_sample_common_names;

    my $id_method = $self->sample_identifier;

    my @all_dna_samples = grep {$_->sample_type =~ /dna/i} $analysis_project->samples;
    my @dna_samples = List::MoreUtils::uniq(@all_dna_samples);
    $self->status_message('Found '. scalar(@dna_samples) .' DNA samples');

    my %dna_samples_by_individual_id;
    for my $dna_sample (@dna_samples) {
        push @{$dna_samples_by_individual_id{$dna_sample->individual->id}}, $dna_sample;
    }
    my @individual_ids = sort keys %dna_samples_by_individual_id;
    $self->status_message('Found '. scalar(@individual_ids) .' indivudals');

    my @subject_mappings;
    for my $individual_id (@individual_ids) {
        my $individual = Genome::Individual->get($individual_id);
        my @individual_dna_samples = @{$dna_samples_by_individual_id{$individual_id}};
        my $individual_dna_sample_count = scalar(@individual_dna_samples);
        $self->status_message('Found '. $individual_dna_sample_count .' samples for individual '. $individual->__display_name__);

        if ( $individual_dna_sample_count > 2) {
            $self->fatal_message('This command currently does not support more than 2 DNA samples per individual.');
        }

        my @normal_samples;
        my @tumor_samples;

        for my $individual_dna_sample (@individual_dna_samples) {
            my $dna_sample_common_name = $individual_dna_sample->common_name || "NULL";
            push(@normal_samples, $individual_dna_sample) if ($dna_sample_common_name =~ /$normal_sample_common_names/i);
            push(@tumor_samples,  $individual_dna_sample) if ($dna_sample_common_name =~ /$tumor_sample_common_names/i);
        }

        if (scalar(@normal_samples) > 1) {
            $self->fatal_message('More than one normal DNA sample was specified for this individual: '. $individual->__display_name__);
        }
        if (scalar(@tumor_samples) > 1) {
            $self->fatal_message('More than one tumor DNA sample was specified for this individual: '. $individual->__display_name__);
        }
        my %data;
        if (@normal_samples && @tumor_samples) {
            %data = (
                tumor_sample => $tumor_samples[0]->$id_method,
                normal_sample => $normal_samples[0]->$id_method,
                snv_variant_list_id => '',
                indel_variant_list_id => '',
                sv_variant_list_id => '',
                tag => 'somatic',
            );
        } elsif (@normal_samples) {
            %data = (
                tumor_sample => $normal_samples[0]->$id_method,
                normal_sample => '',
                snv_variant_list_id => '',
                indel_variant_list_id => '',
                sv_variant_list_id => '',
                tag => 'germline',
            );
        } elsif (@tumor_samples) {
            %data = (
                tumor_sample => $tumor_samples[0]->$id_method,
                normal_sample => '',
                snv_variant_list_id => '',
                indel_variant_list_id => '',
                sv_variant_list_id => '',
                tag => 'tumor-only',
            );
        }
        # TODO: See if subject mapping already exists before adding to the list
        push @subject_mappings, \%data;
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->file_path,
        separator => "\t",
        headers => \@headers,
        print_headers => 0,
    );

    for my $subject_mapping (@subject_mappings) {
        $writer->write_one($subject_mapping);
    }
    $writer->output->close();

    my $count = scalar(@subject_mappings);
    $self->status_message('Predicted '. $count .' new subject mappings.');

    return $count;
}



1;
