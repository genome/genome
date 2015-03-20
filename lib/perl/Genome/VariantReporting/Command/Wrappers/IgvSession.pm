package Genome::VariantReporting::Command::Wrappers::IgvSession;

use strict;
use warnings;
use Genome;
use URI;

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
            is => 'Genome::VariantReporting::Framework::Component::Report::MergedReport',
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

    my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlBasic->create(
        output_file    => $self->_temp_file_path,
        resource_files => $self->resource_files,
        reference_name => $self->igv_reference_name,
    );
    unless ($dumpXML->execute) {
        confess $self->error_message("Failed to create IGV xml file");
    }
}

sub bams {
    my $self = shift;
    return $_JSON_CODEC->decode($self->bam_hash_json);
}

sub uri {
    my $file = shift;
    return URI->new_abs($file, $ENV{GENOME_SYS_SERVICES_FILES_URL})->as_string;
}

sub resource_files {
    my $self = shift;

    my %bams = %{$self->bams};
    my %reference_files = map { uri($bams{$_}) => $_ } keys %bams;


    for my $bed_report ($self->merged_bed_reports) {
        my @report_users = map { $_->users('label like' => 'report:%') } $bed_report->report_results;
        my ($process) = grep { $_->isa('Genome::VariantReporting::Process::Trio') } $bed_report->children;
        $reference_files{uri($bed_report->report_path)} = bed_file_label($bed_report->category, $process) || $bed_report->report_path;
    }

    return \%reference_files;
}

sub bed_file_label {
    my ($category, $process) = @_;

    return unless defined($process);

    my %labels = (
        docm      => 'Recurrent AML Variants',
        followup  => sprintf('Followup(%s) Variants', $process->followup_sample->name),
        discovery => sprintf('Discovery(%s) Variants', $process->tumor_sample->name),
        germline  => sprintf('Germline(%s) Variants', $process->normal_sample->name)
    );

    return $labels{$category};
}

sub igv_reference_name {
    my $self = shift;
    my $reference_sequence_name_cmd = Genome::Model::Tools::Analysis::ResolveIgvReferenceName->execute(
        reference_name => $self->reference_name,
    );
    return $reference_sequence_name_cmd->igv_reference_name;
}
1;

