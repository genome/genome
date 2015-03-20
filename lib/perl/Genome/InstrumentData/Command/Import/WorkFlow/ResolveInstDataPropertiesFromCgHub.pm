package Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataPropertiesFromCgHub;

use strict;
use warnings;

use Genome;

use Genome::Model::Tools::CgHub::Query;
require File::Spec;
require File::Temp;
use Path::Class;

class Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataPropertiesFromCgHub { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::ResolveInstDataProperties',
    has_optional_transient => {
        metadata_file => { is => 'Text', },
    },
};

sub execute {
    my $self = shift;

    my $resolve_instdata_props = $self->_resolve_instrument_data_properties;
    return if not $resolve_instdata_props;

    my $metadata_file = $self->_resolve_metadata_file;
    return if not $metadata_file;

    my $add_properties_from_metadata = $self->_add_properties_from_metadata;
    return if not $add_properties_from_metadata;

    return 1;
}

sub _resolve_metadata_file {
    my $self = shift;
    $self->debug_message('Resolve metadata file...');

    # check if DL'd
    my $source = $self->source;
    my $source_dir = Path::Class::file($source)->dir;
    my $metadata_file_basename = 'metadata.xml';
    my $metadata_file = File::Spec->join($source_dir, $metadata_file_basename);
    return $self->metadata_file($metadata_file) if -s $metadata_file;

    # DL to tmpdir
    # Need UUID to DL
    my $uuid = delete $self->resolved_instrument_data_properties->{uuid};
    if ( not defined $uuid ) {
        $self->warning_message('No uuid found in instrument data properties! It is required to download the metadata file.');
        return 1; # ok for now
    }

    # to tmpdir
    my $tempdir = File::Temp::tempdir(CLEANUP => 1);
    if ( not $tempdir ) {
        $self->error_message('Failed to create tmp dir!');
        return;
    }
    $metadata_file = File::Spec->join($tempdir, $metadata_file_basename);

    my $query = Genome::Model::Tools::CgHub::Query->execute(
        uuid => $uuid,
        xml_file => $metadata_file,
    );
    if ( not $query->result ) {
        $self->error_message('Failed to execute cg hub query!');
        return;
    }
    $self->metadata_file($metadata_file);
    $self->debug_message("Metadata file: $metadata_file");

    $self->debug_message('Resolve metadata file...done');
    return 1;
}

sub _add_properties_from_metadata {
    my $self = shift;;
    $self->debug_message('Add properties from metadata...');

    my $metadata_file = $self->metadata_file;
    return 1 if not -s $metadata_file;

    my $metadata = Genome::Model::Tools::CgHub::Metadata->create(
        metadata_file => $metadata_file,
    );

    my %args = (
        aliquot_id => 0,
        analysis_id => 0,
        participant_id => 0,
        sample_id => 0,
        import_source_name => 1,
        tcga_name => 1,
    );
    my $properties = $self->resolved_instrument_data_properties;
    for my $arg_name (keys %args){
        if(defined $properties->{$arg_name}){
            $self->status_message('Argument (%s) was passed in as (%s)', $arg_name, $properties->{$arg_name});
            next;
        }
        my $meta_value = ( $metadata ? $metadata->get_attribute_value($arg_name) : undef );
        if ( defined $meta_value ) {
            $self->status_message('Argument (%s) was found in the metadata file in as (%s)', $arg_name, $meta_value);
            $properties->{$arg_name} = $meta_value;
            next;
        }
        if ( $args{$arg_name} ) { # required
            die $self->error_message("Required argument ($arg_name) was not found in metadata file and was not passed in via instrument data properties.");
        }
    }

    $properties->{description} ||= "TCGA Bam imported from CG Hub for ".$properties->{tcga_name};

    if ( not $properties->{target_region_set_name} and $metadata ) {
        $properties->{target_region_set_name} = $metadata->target_region;
    }

    if ( not $properties->{target_region_set_name} ) {
        die $self->error_message('Target region is requried!');
    }

    unless ($properties->{target_region_set_name} eq 'none') {
        my @feature_lists = Genome::FeatureList->get(name => $properties->{target_region_set_name});
        if (not @feature_lists or @feature_lists > 1) {
            $self->error_message("Invalid target region: " . $properties->{target_region_set_name});
            return;
        }
    }

    # If these do not exist that may be ok.
    my $bam_file_size = $self->_bam_file_size_from_metadata($metadata);
    $properties->{source_size} = $bam_file_size if $bam_file_size;

    my $source_md5 = $self->_bam_file_md5_from_metadata($metadata);
    $properties->{source_md5} = $source_md5 if $source_md5;

    $self->debug_message('Add properties from metadata...');
    return $self->resolved_instrument_data_properties($properties);
}

sub _bam_file_size_from_metadata {
    my ($self, $metadata) = @_;

    Carp::confess('No metadata given to get bam file size!') if not $metadata;

    my $bam_file_name = $metadata->bam_file_names;
    return if not $bam_file_name;

    return $metadata->file_size_for_file_name($bam_file_name);
}

sub _bam_file_md5_from_metadata { 
    my ($self, $metadata) = @_;

    Carp::confess('No metadata given to get bam file size!') if not $metadata;

    my $bam_file_name = $metadata->bam_file_names;
    return if not $bam_file_name;

    my $checksum_type = $metadata->checksum_type_for_file_name($bam_file_name);
    return if not $checksum_type or lc($checksum_type) ne 'md5';

    return $metadata->checksum_content_for_file_name($bam_file_name);
}

1;

