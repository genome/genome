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
        genotype_build_id => {
            is => 'Text',
        },
        on_target_list => {
            is => "Genome::FeatureList",
            is_optional => 1,
        },
    ],
    has_param => [
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

    return $self;
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

    $self->_error("Could not find instrument data for id ".$self->aligned_bam_result_id) unless $bam_result;
    my $path = $bam_result->bam_path;
    unless (-s $path) {
        $self->_error("Could not get bam file for ".$bam_result->id);
    }
    return $path;
}

sub _resolve_vcf_file {
    my $self = shift;
    my $genotype_build = Genome::Model::Build::GenotypeMicroarray->get(id => $self->genotype_build_id);

    $self->_error("Could not find genotype build for id ".$self->genotype_build_id) unless $genotype_build;
    my $vcf = Genome::Sys->create_temp_file_path;
    my $params = {
        build => $genotype_build,
        variation_list_build => $genotype_build->dbsnp_build,
        output => $vcf,
    };
    if ($genotype_build->reference_sequence_build->allosome_names) {
        $params->{filters} = ["chromosome:exclude=".$genotype_build->reference_sequence_build->allosome_names];
    }
    my $rv = Genome::Model::GenotypeMicroarray::Command::ExtractToVcf->execute($params);
    unless ($rv and -s $vcf) {
        $self->_error("Could not get vcf file for ".Data::Dumper::Dumper($params));
    }
    return $self->_clean_vcf($vcf);
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

