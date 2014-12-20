package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub;

use strict;
use warnings;

use Genome;

require Genome::Model::Tools::CgHub::GeneTorrent;
require Genome::Model::Tools::CgHub::Metadata;
require Genome::Model::Tools::CgHub::Query;
require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
    has => {
        lsf_resource => {
            default_value => Genome::Model::Tools::CgHub::GeneTorrent->__meta__->property_meta_for_name('lsf_resource')->default_value,
        },
    },
    has_calculated => {
        metadata_file => {
            calculate_from => [qw/ working_directory /],
            calculate => q| return File::Spec->join($working_directory, 'metadata.xml'); |,
        },
        source_path_basename => {
            calculate_from => [qw/ source_path /],
            calculate => q| return File::Basename::basename($source_path).'.bam'; |,
        },
        uuid => {
            calculate_from => [qw/ source_path /],
            calculate => q| return File::Basename::basename($source_path); |,
        },
    },
};

sub _source_path_size {
    my $self = shift;

    my $metadata = $self->_metadata;
    return if not $metadata;

    my $bam_file_name = $metadata->bam_file_names;
    if ( not defined $bam_file_name ) {
        $self->error_message('No source path found in CG Hub metadata!');
        return;
    }

    my $source_path_sz = $metadata->file_size_for_file_name($bam_file_name);

    if ( not defined $source_path_sz ) {
        $self->error_message("No source path size found in CG Hub metadata for $bam_file_name!");
        return;
    }

    return $source_path_sz;
}

sub _retrieve_source_path {
    my $self = shift;

    my $gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->execute(
        uuid => $self->uuid,
        target_path => $self->destination_path,
    );
    if ( not $gene_torrent->result ) {
        $self->error_message('Failed to execute cg hub gene torrent!');
        return;
    }

    return 1;
}

sub source_md5 {
    my $self = shift;

    my $metadata = $self->_metadata;
    return if not $metadata;

    my $bam_file_name = $metadata->bam_file_names;
    my $checksum_type = $metadata->checksum_type_for_file_name($bam_file_name);
    if ( not $checksum_type or lc($checksum_type) ne 'md5' ) {
        $self->debug_message('Checksum from metadata is not MD5, skipping storing it.');
        return;
    }

    return $metadata->checksum_content_for_file_name($bam_file_name);
}

sub _metadata {
    my $self = shift;

    return $self->{_metadata} if $self->{_metadata};

    my $retrieve_ok = $self->_retrieve_metadata_path;
    return if not $retrieve_ok;

    $self->{_metadata} = Genome::Model::Tools::CgHub::Metadata->create(
        metadata_file => $self->metadata_file,
    );
    if ( not $self->{_metadata} ) {
        $self->error_message('Failed to load metadata from file! '.$self->metadata_file);
        return;
    }

    return $self->{_metadata};
}

sub _retrieve_metadata_path {
    my $self = shift;
    $self->debug_message('Retrieve metadata path from CG Hub...');

    my $uuid = $self->uuid;
    $self->debug_message("UUID: $uuid");
    my $metadata_file = $self->metadata_file;
    $self->debug_message("Metadata file: $metadata_file");

    my $query = Genome::Model::Tools::CgHub::Query->execute(
        uuid => $uuid,
        xml_file => $metadata_file,
    );
    if ( not $query->result ) {
        $self->error_message('Failed to execute cg hub query!');
        return;
    }

    if ( not -s $metadata_file ) {
        $self->error_message("Successfully executed query, but the metadata path does not exist!");
        return;
    }

    $self->debug_message('Retrieve metadata path from CG Hub...done');
    return 1;
}

1;

