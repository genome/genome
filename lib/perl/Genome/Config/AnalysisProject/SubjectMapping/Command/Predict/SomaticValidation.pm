package Genome::Config::AnalysisProject::SubjectMapping::Command::Predict::SomaticValidation;

use strict;
use warnings;

use Genome;
use Memoize;
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
    has_transient_optional => [
        existing_subject_mapping_set => {
            is => 'Set::Scalar',
            doc => 'The set of existing subject mappings.',
        },
        new_subject_mapping_set => {
            is => 'Set::Scalar',
            doc => 'The set of new subject mappings.',
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

    my $samples_by_individual_id = $self->resolve_samples_by_individual_id;
    my @individual_ids = sort keys %{$samples_by_individual_id};
    $self->status_message('Found '. scalar(@individual_ids) .' individuals');

    for my $individual_id (@individual_ids) {
        my $individual = Genome::Individual->get($individual_id);

        my @samples = @{$samples_by_individual_id->{$individual_id}};

        my $sample_set = Set::Scalar->new(@samples);
        $self->status_message('Found '. $sample_set->size .' samples for individual '. $individual->__display_name__);

        my @normal_samples;
        my @tumor_samples;

        for my $sample ($sample_set->members) {
            my $sample_common_name = $sample->common_name || 'NULL';
            push(@normal_samples, $sample) if ($sample_common_name =~ /$normal_sample_common_names/i);
            push(@tumor_samples,  $sample) if ($sample_common_name =~ /$tumor_sample_common_names/i);
        }

        my $normal_sample_set = Set::Scalar->new(@normal_samples);
        my $tumor_sample_set = Set::Scalar->new(@tumor_samples);

        if ($normal_sample_set->size > 1) {
            $self->fatal_message('The following samples for individual '. $individual->__display_name__ .' are ALL found to be normal: '. join(',', map {$_->$id_method} $normal_sample_set->members));
        }

        my $unique_sample_set = $sample_set->unique($normal_sample_set,$tumor_sample_set);
        if ($unique_sample_set) {
            $self->fatal_message('The following samples for individual: '. $individual->__display_name__ .' do not match either tumor or normal criteria: '. join(',', map{$_->$id_method} $unique_sample_set->members));
        }

        my $intersection_sample_set = $normal_sample_set->intersection($tumor_sample_set);
        if ($intersection_sample_set) {
            $self->fatal_message('The following samples for individual '. $individual->__display_name__ .' meet both tumor and normal match criteria: '. join(',', map {$_->$id_method} $intersection_sample_set->members));
        }

        if ($normal_sample_set && $tumor_sample_set) {
            my ($normal_sample) = $normal_sample_set->members();
            for my $tumor_sample ($tumor_sample_set->members) {
                $self->add_subject_mapping($tumor_sample,$normal_sample,undef,undef,undef,'discovery');
            }
        } elsif ($normal_sample_set) {
            my ($normal_sample) = $normal_sample_set->members();
            $self->add_subject_mapping($normal_sample,undef,undef,undef,undef,'germline');
        } elsif ($tumor_sample_set) {
            for my $tumor_sample ($tumor_sample_set->members) {
                $self->add_subject_mapping($tumor_sample,undef,undef,undef,undef,'tumor-only');
            }
        }
    }
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->file_path,
        separator => "\t",
        headers => \@headers,
        print_headers => 0,
    );

    for my $subject_mapping ($self->new_subject_mapping_set->members) {
        $writer->write_one($subject_mapping);
    }
    $writer->output->close();

    $self->status_message('Predicted '. $self->new_subject_mapping_set->size .' new subject mappings');

    return 1;
}


sub existing_subject_mapping_set {
    my $self = shift;

    unless (defined($self->__existing_subject_mapping_set)) {

        my $id_method = $self->sample_identifier;
        my @subject_mappings;
        for my $subject_mapping ($self->analysis_project->subject_mappings) {
            my ($tumor_sample, $normal_sample, $snv_id, $indel_id, $sv_id, $tag_name);
            for my $subject_bridge ($subject_mapping->subject_bridges) {
                if ($subject_bridge->label eq 'tumor_sample') {
                    $tumor_sample = $subject_bridge->subject;
                }
                if ($subject_bridge->label eq 'normal_sample') {
                    $normal_sample = $subject_bridge->subject;
                }
            }
            for my $input ($subject_mapping->inputs) {
                if ($input->key eq 'snv_variant_list_id') {
                    $snv_id = $input->value;
                }
                if ($input->key eq 'indel_variant_list_id') {
                    $indel_id = $input->value;
                }
                if ($input->key eq 'sv_variant_list_id') {
                    $sv_id = $input->value;
                }
            }
            for my $tag ($subject_mapping->tags) {
                if (defined($tag_name)) {
                    $self->fatal_message('Unable to handle multiple subject mapping tags!');
                }
                $tag_name = $tag->name;
            }
            my $data = $self->_subject_mapping_data_hash_ref($tumor_sample,$normal_sample,$snv_id,$indel_id,$sv_id,$tag_name);
            push @subject_mappings, $data;
        }

        my $subject_mapping_set = Set::Scalar->new(@subject_mappings);
        $self->status_message('Found '. $subject_mapping_set->size .' existing subject mappings.');
        $self->__existing_subject_mapping_set($subject_mapping_set);
    }
    return $self->__existing_subject_mapping_set;
}

sub new_subject_mapping_set {
    my $self = shift;
    unless (defined($self->__new_subject_mapping_set)) {
        my $new_subject_mapping_set = Set::Scalar->new();
        $self->__new_subject_mapping_set($new_subject_mapping_set);
    }
    return $self->__new_subject_mapping_set;
}

sub add_subject_mapping {
    my $self = shift;

    my $data = $self->_subject_mapping_data_hash_ref(@_);

    if ($self->existing_subject_mapping_set->element($data)) {
        $self->warning_message('Skipping existing subject mapping betweern tumor \''. $data->{tumor_sample} .'\' and normal \''. $data->{normal_sample} .'\'!');
        return;
    }
    if ($self->new_subject_mapping_set->element($data)) {
        $self->warning_message('New subject mapping already exists in set, skipping!');
        return;
    }
    $self->new_subject_mapping_set->insert($data);
    return 1;
}

sub resolve_samples_by_individual_id {
    my $self = shift;

    my @all_dna_samples = grep {$_->sample_type =~ /dna/i && !($_->is_rna) } $self->analysis_project->samples;
    my @dna_samples = List::MoreUtils::uniq(@all_dna_samples);
    $self->status_message('Found '. scalar(@dna_samples) .' DNA samples');

    my %samples_by_individual_id;
    for my $dna_sample (@dna_samples) {
        push @{$samples_by_individual_id{$dna_sample->individual->id}}, $dna_sample;
    }
    return \%samples_by_individual_id
}

memoize '_subject_mapping_data_hash_ref';

sub _subject_mapping_data_hash_ref {
    my $self = shift;
    my ($tumor, $normal, $snv_id, $indel_id, $sv_id, $tag_name) = @_;
    my $id_method = $self->sample_identifier;
    my %data = (
        tumor_sample  => $tumor->$id_method,
        normal_sample => $normal ? $normal->$id_method : '',
        snv_variant_list_id => $snv_id || '',
        indel_variant_list_id => $indel_id || '',
        sv_variant_list_id => $sv_id || '',
        tag => $tag_name,
    );
    return \%data;
}

1;
