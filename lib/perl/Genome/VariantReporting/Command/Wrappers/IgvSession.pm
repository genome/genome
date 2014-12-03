package Genome::VariantReporting::Command::Wrappers::IgvSession;

use strict;
use warnings;
use Genome;

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::VariantReporting::Command::Wrappers::IgvSession {
    is => 'Genome::SoftwareResult::StageableSimple',
    has_input => [
        bam_hash_json_lookup => {
            is => 'String',
            doc => 'Hash where keys are the labels for the bam track and the values are the bam file paths',
        },
        genome_name => {
            is => 'String',
            doc => 'Name for the file and tracks',
        },
        merged_bed_reports => {
            is => 'Genome::VariantReporting::Framework::MergedReport',
            is_many => 1,
            doc => 'Reports to get the bed files from',
        },
        reference_name => {
            is => 'String',
            doc => 'name of the reference sequence build',
        },
    ],
    has_transient_optional => [
        bam_hash_json => {
            is => 'String',
            doc => 'Hash where keys are the labels for the bam track and the values are the bam file paths',
        },
    ],
};

sub output_file_name {
    return "igv.xml";
}

sub output_file_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->output_file_name);
}

sub _temp_file_path {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_file_name);
}

sub _run {
    my $self = shift;

    my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
        bams             => $self->bam_paths,
        labels           => $self->bam_labels,
        output_file      => $self->_temp_file_path,
        genome_name      => $self->genome_name,
        review_bed_files => $self->bed_files,
        reference_name   => $self->igv_reference_name,
    );
    unless ($dumpXML->execute) {
        confess $self->error_message("Failed to create IGV xml file");
    }
}

sub bams {
    my $self = shift;
    return $_JSON_CODEC->decode($self->bam_hash_json);
}

sub bam_paths {
    my $self = shift;
    #return join(',', map {File::Spec->join($ENV{GENOME_SYS_SERVICES_FILES_URL}, $_)} values %{$self->bams});
    return join(',', values %{$self->bams});
}

sub bam_labels {
    my $self = shift;
    return join(',', keys %{$self->bams});
}

sub bed_files {
    my $self = shift;
    return [map {$_->report_path} $self->merged_bed_reports];
}

sub igv_reference_name {
    my $self = shift;
    my $reference_sequence_name_cmd = Genome::Model::Tools::Analysis::ResolveIgvReferenceName->execute(
        reference_name => $self->reference_name,
    );
    return $reference_sequence_name_cmd->igv_reference_name;
}
1;

