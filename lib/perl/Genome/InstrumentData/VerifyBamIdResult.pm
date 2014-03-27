package Genome::InstrumentData::VerifyBamIdResult;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

class Genome::InstrumentData::VerifyBamIdResult {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        aligned_bam_result_id => {
            is => 'Text',
        },
        on_target_list => {
            is => "Genome::FeatureList",
            is_optional => 1,
        },
        sample => {
            is => "Genome::Sample",
        },
        known_sites_build => {
            is => "Genome::Model::Build::ImportedVariationList",
        },
    ],
    has_param => [
        genotype_filters => {
            is => 'Text',
            is_many => 1,
        },
        max_depth => {
            is => "Integer",
        },
        precise => {
            is => 'Boolean',
        },
        version => {
            is => "Text",
        },
    ],
    has_metric => [
        freemix => {
            is => "UR::Value::Number",
        },
        chipmix => {
            is => "UR::Value::Number",
        },
    ],
};

sub _error {
    my ($self, $msg) = @_;
    $self->error_message($msg);
    $self->delete;
    die $self->error_message;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return if not $self;
    $self->_error("Failed to prepare staging directory") unless $self->_prepare_staging_directory;

    $self->_error("Failed to run verifyBamID") unless $self->_run_verify_bam_id;


    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;
    
    $self->_error("Failed to add metrics") unless $self->_add_metrics;

    return $self;
}

sub _add_metrics {
    my $self = shift;
    my $self_sm = File::Spec->join($self->output_dir, "output.selfSM");
    my $in = Genome::Utility::IO::SeparatedValueReader->create(
        separator => "\t",
        input => $self_sm,
    );
    my $metrics = $in->next;
    $self->freemix($metrics->{FREEMIX});
    $self->chipmix($metrics->{CHIPMIX});
    return 1;
}

sub _run_verify_bam_id {
    my $self = shift;

    my $bam_file = $self->_resolve_bam_file;
    my $vcf_file = $self->_resolve_vcf_file;
    my $out_prefix = File::Spec->join($self->temp_staging_directory, "output");
    return Genome::Model::Tools::VerifyBamId->execute(vcf => $vcf_file,
                bam => $bam_file, out_prefix => $out_prefix, max_depth => $self->max_depth,
                precise => $self->precise, version => $self->version);
}

sub _resolve_bam_file {
    my $self = shift;
    my $bam_result = Genome::InstrumentData::AlignedBamResult->get($self->aligned_bam_result_id);

    $self->_error("Could not find alignment result for id ".$self->aligned_bam_result_id) unless $bam_result;
    my $path = $bam_result->bam_path;
    unless (-s $path) {
        $self->_error("Could not get bam file for ".$bam_result->id);
    }
    return $path;
}

sub _resolve_vcf_file {
    my $self = shift;
    my $genotype_vcf_result = $self->_resolve_genotype_vcf_result;
    my $vcf = $genotype_vcf_result->vcf_path;
    unless (-s $vcf) {
        $self->_error("Could not get vcf file for genotype vcf".$genotype_vcf_result);
    }
    return $self->_clean_vcf($vcf);
}

sub _resolve_genotype_vcf_result {
    my $self = shift;

    my %params = (
        sample => $self->sample,
        known_sites_build => $self->known_sites_build,
    );
    if ($self->genotype_filters) {
        $params{filters} = [$self->genotype_filters];
    }
    my $result = Genome::InstrumentData::GenotypeVcf->get_or_create(%params);
    $self->_error("Could not get or create genotype vcf result") unless $result;
    return $result;
}

sub _clean_vcf {
    my $self = shift;
    my $vcf_path = shift;

    if ($self->on_target_list) {
        $self->debug_message("Using on_target_list");
        my $on_target_path = Genome::Sys->create_temp_file_path;
        my $on_target_bed = $self->on_target_list->processed_bed_file;
        my $rv = Genome::Model::Tools::BedTools::Intersect->execute(
            input_file_a => $vcf_path,
            input_file_b => $on_target_bed,
            input_file_a_format => "bed",
            intersection_type => "unique",
            output_file => $on_target_path,
            header => 1,
        );
        $self->_error("Could not intersect with on target bed") unless $rv;
        $vcf_path = $on_target_path;
    }

    $self->debug_message("Using cleaned vcf at $vcf_path");
    return $vcf_path;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("verifybamidresult-%s-%s-%s-%s",           $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}
1;

