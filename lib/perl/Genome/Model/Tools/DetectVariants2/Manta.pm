package Genome::Model::Tools::DetectVariants2::Manta;

use warnings;
use strict;

use Genome;

use File::Spec qw(join);
use Data::Dumper;

my $CPU = 12;
my $RAM = 8;

class Genome::Model::Tools::DetectVariants2::Manta {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
    has_param => [
        lsf_resource => {
            default_value => Genome::Config::get('lsf_resource_dv2_manta_sv'),
        },
    ],
};

sub _cpu {
    return $CPU;
}

sub _ram {
    return $RAM;
}

sub _detect_variants {
    my $self = shift;

    my $tumor_ar = $self->aligned_reads_input;
    my $normal_ar = $self->control_aligned_reads_input;

    my $tumor_bam_file = $self->aligned_reads_input->bam_file;
    my $normal_bam_file = $self->control_aligned_reads_input->bam_file;

    my $working_directory = $self->_temp_staging_directory;

    my %config_params = (
        tumor_bam_file => $tumor_bam_file,
        normal_bam_file => $normal_bam_file,
        version => $self->version,
        working_directory => $working_directory,
        reference_fasta => $self->reference_sequence_input,
    );

    # TODO : Determine if any of the input instrument data is exome, rna (if rna, unstranded?)

    my $config = Genome::Model::Tools::Manta::Config->create(%config_params);

    unless ($config) {
        $self->fatal_message('Failed to create Manta config command with params: '. Data::Dumper::Dumper(%config_params) );
    }
    unless ($config->execute()) {
        $self->fatal_message('Could not execute Manta config command!');
    }

    my %run_params = (
        working_directory => $working_directory,
        jobs => $self->_cpu,
        memory => $self=>_ram,
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
