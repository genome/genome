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
            default_value => Genome::Config::get('lsf_resource_dv2_manta'),
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

    my $working_directory = $self->_temp_staging_directory;

    my %config_params = (
        tumor_bam_file => $self->aligned_reads_input,
        normal_bam_file => $self->control_aligned_reads_input,
        version => $self->version,
        working_directory => $working_directory,
        reference_fasta => $self->reference_sequence_input,
    );

    # TODO : Determine if this is exome or rna (if rna, unstranded?)
    # One option is to define these as detector params in the processing profile, ie. $self->params
    # If so, parse the param string to make the necessary param to hash ref translation
    # The same could be done with the config_file option but it will have a value rather than a simple flag
    # Example : --exome --config-file=/my/foo/config.txt
    # ( exome => 1, config_file => '/my/foo/config/txt' )

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
        memory => $self->_ram,
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
