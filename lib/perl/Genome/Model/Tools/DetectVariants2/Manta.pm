package Genome::Model::Tools::DetectVariants2::Manta;

use warnings;
use strict;

use Genome;

use Genome::Sys::LSF::ResourceParser qw(parse_lsf_params);

use File::Spec qw(join);
use Text::ParseWords qw(shellwords);
use Data::Dumper;

class Genome::Model::Tools::DetectVariants2::Manta {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
    has_param => [
        lsf_resource => {
            default_value => Genome::Config::get('lsf_resource_dv2_manta'),
        },
    ],
};

sub _gb_ram {
    my $self = shift;

    my $lsf_params = parse_lsf_params($self->lsf_resource);

    my $kb_ram = $lsf_params->{'rLimits'}->{'RSS'};

    if (!defined($kb_ram)) {
        $self->fatal_message('Unable to resolve RAM from the lsf resource string: '. $self->lsf_resource);
    }

    return ($kb_ram / 1_000_000);
}

sub _detect_variants {
    my $self = shift;

    my $working_directory = $self->_temp_staging_directory;

    # Optional params defined as detector params in the processing profile
    # parse the param string to make the necessary param to hash ref translation
    # Example : --exome --config-file=/my/foo/config.txt
    # ( exome => 1, config_file => '/my/foo/config/txt' )
    # Also see Genome::Model::Tools::DetectVariants2::Mutect::parse_params
    my @resolved_config_params = (
        '--tumor-bam-file' => $self->aligned_reads_input,
        '--version' => $self->version,
        '--working-directory' => $working_directory,
        '--reference-fasta' => $self->reference_sequence_input,
    );

    # Determine which manta output VCF should be used as the "high-quality" result
    # If single-sample, assume tumor-only (NOT GERMLINE).  We know somatic if a normal is defined.
    my $manta_vcf_filename = 'tumorSV.vcf.gz';

    if ( defined($self->control_aligned_reads_input) ) {
        push @resolved_config_params, ('--normal-bam-file' => $self->control_aligned_reads_input);
        $manta_vcf_filename = 'somaticSV.vcf.gz';
    }

    my @dv2_config_params = Text::ParseWords::shellwords($self->params);
    my ($config_cmd_class, $config_cmd_params) = Genome::Model::Tools::Manta::Config->resolve_class_and_params_for_argv(@dv2_config_params,@resolved_config_params);
    my $config_cmd = Genome::Model::Tools::Manta::Config->create(%$config_cmd_params);

    unless ($config_cmd) {
        $self->fatal_message('Failed to create Manta config command with params: '. Data::Dumper::Dumper($config_cmd_params) );
    }
    unless ($config_cmd->execute()) {
        $self->fatal_message('Could not execute Manta config command!');
    }

    my %run_params = (
        working_directory => $working_directory,
        jobs => Genome::SoftwareResult->_available_cpu_count,
        memory => $self->_gb_ram,
    );

    my $run = Genome::Model::Tools::Manta::Run->create(%run_params);

    unless ($run) {
        $self->fatal_message('Failed to create Manta run command with params: '. Data::Dumper::Dumper(%run_params));
    }
    unless ($run->execute) {
        $self->fatal_message('Could not execute Manta run command!');
    }

    # Move the VCF to where the DV2 Dispatcher and the DetectVariants build command expect output with a basename like 'svs.hq'
    my $manta_vcf_path = File::Spec->join($working_directory,'results/variants/'. $manta_vcf_filename);
    if (-e $manta_vcf_path) {
        my $sv_staging_path = $self->_sv_staging_output .'.vcf.gz';
        # Could use move_file if disk where a concern, but retaining the original results/ directory contents from manta
        Genome::Sys->copy_file($manta_vcf_path, $sv_staging_path);
        Genome::Sys->copy_file($manta_vcf_path .'.tbi', $sv_staging_path .'.tbi');
    } else {
        $self->fatal_message('Failed to find the Manta somatic SV VCF file expected at: '. $manta_vcf_path);
    }
}

sub has_version {
    my $self    = shift;
    my $version = shift;

    my $path = Genome::Model::Tools::Manta::Base->path_for_version($version);
    if (-e $path) {
        return 1;
    }
    else {
        return 0;
    }
}

sub _sort_detector_output {
    # This command does not produce a BED output file
    return 1;
}


1;
