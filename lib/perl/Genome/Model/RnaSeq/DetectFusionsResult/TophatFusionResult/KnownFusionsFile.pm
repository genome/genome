package Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult::KnownFusionsFile;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::TophatFusionResult::KnownFusionsFile {
    is => 'Genome::SoftwareResult::Stageable',
    has => [
        fusion_file_path => {
            is => "Text",
            is_metric => 1,
            doc => 'path on network disk to the file to be used'
        },
        file_source_description => {
            is => "Text",
            is_input => 1,
            doc => "brief description of where this file came from or how it was made"
        },
        mcl_file => {
            is => "Text",
            doc => "the mcl file for this known fusions result",
            is_calculated => 1,
            calculate => '$output_dir."/mcl"',
            calculate_from => ["output_dir"],
        }
    ],
    doc => 'file used by tophat fusion for known fusion annotation'
};


sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/known-fusions-file/' . $self->id;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_) or return;
    $self->_prepare_staging_directory;

    Genome::Sys->copy_file($self->fusion_file_path, $self->temp_staging_directory . "/mcl");

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub __display_name__ {
    my $self = shift;
    return $self->file_source_description . " (" . $self->id . ")";
}

1;
