package Genome::Model::Tools::DetectVariants2::PindelRegion;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::PindelRegion {
    is => 'Genome::Model::Tools::DetectVariants2::Detector',
};


sub _detect_variants {
    my $self = shift;

    my ($version, $param_str) = map{$self->$_}qw(version params);
    my @compatible_versions = Genome::Model::Tools::Pindel::RunPindel->pindel_region_compatible_versions;
    my $compatible_versions = join ',', @compatible_versions;

    unless (grep{$version eq $_}@compatible_versions) {
        $self->fatal_message("pindel version $version does not have option to run region file. Compatible versions are $compatible_versions");
    }

    my $region_file;

    if (-s $param_str) {
        $region_file = $param_str;
    }
    else {
        my $feature_list = Genome::FeatureList->get($param_str);
        if ($feature_list) {
            $region_file = $feature_list->file_path;
        }
        unless ($region_file and -s $region_file) {
            $self->fatal_message("Fail to get region file from $param_str");
        }
    }

    my %params = (
        aligned_reads_input => $self->aligned_reads_input,
        reference_build_id  => $self->reference_build_id,
        output_directory    => $self->_temp_staging_directory,
        version             => $version,
        region_file         => $region_file,
    );

    $params{control_aligned_reads_input} = $self->control_aligned_reads_input if $self->control_aligned_reads_input;

    my $cmd = Genome::Model::Tools::Pindel::RunPindel->create(%params);
    
    unless ($cmd) {
        $self->fatal_message('Failed to create pindel command');
    }
    unless($cmd->execute()) {
        $self->fatal_message('Failed to execute pindel');
    }
    
    return 1;
}


sub has_version {
    my ($self, $version) = @_;

    $version = $self->version unless defined $version;
    my @versions = Genome::Model::Tools::Pindel::RunPindel->pindel_region_compatible_versions;

    return 1 if grep{$version eq $_}@versions;
    return 0;
}
  

sub _sort_detector_output {
    # This command does not produce a BED output file
    return 1;
}

1;
