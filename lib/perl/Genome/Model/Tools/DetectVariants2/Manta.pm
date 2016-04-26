package Genome::Model::Tools::DetectVariants2::Manta;

use warnings;
use strict;

use Genome;

use Genome::Sys::LSF::ResourceParser qw(parse_lsf_params);

use File::Spec qw(join);
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
        '--normal-bam-file' => $self->control_aligned_reads_input,
        '--version' => $self->version,
        '--working-directory' => $working_directory,
        '--reference-fasta' => $self->reference_sequence_input,
    );

    my @dv2_config_params = split(' ',$self->params);
    my ($config_cmd_class, $config_cmd_params) = Genome::Model::Tools::Manta::Config->resolve_class_and_params_for_argv(@dv2_config_params,@resolved_config_params);

    my $config_cmd = Genome::Model::Tools::Manta::Config->create($config_cmd_params);

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

    # Make a symlink to the actual VCF as the DV2 Dispatcher expects an output of svs.hq rather than BED or VCF output
    my $vcf_path = File::Spec->join($working_directory,'results/variants/somaticSV.vcf.gz');
    if (-e $vcf_path) {
        symlink($vcf_path,$self->_sv_staging_output);
    } else {
        $self->fatal_message('Failed to find the Manta somatic SV VCF file expected at: '. $vcf_path);
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


1;
