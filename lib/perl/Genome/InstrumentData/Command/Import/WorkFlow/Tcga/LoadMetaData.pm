package Genome::InstrumentData::Command::Import::WorkFlow::Tcga::LoadMetaData;

use strict;
use warnings;

use Genome;

use XML::Simple;

class Genome::InstrumentData::Command::Import::WorkFlow::Tcga::LoadMetaData { 
    is => 'Command::V2',
    has_input => {
        metadata_path => {
            is => "Text",
            is_optional => 1,
            doc => "The fullpath of a .xml metadata file which contains (among other things) the md5 of the bam file.  If not specified, a metadata file is looked for in the same directory as the bam file (named metadata.xml).",
        },
    },
    has_output => {
        instrument_data_attributes => {
            is => 'Hash',
            default_value => {},
            doc => 'Instrument data attributes retrieved from the meta data XML.',
        },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my $metadata_path = $self->metadata_path;
    my $metadata_path_ok = eval{ Genome::Sys->validate_file_for_reading($metadata_path); };
    if ( not $metadata_path_ok ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ metadata_path /],
            desc => 'Metadata path is not a valid file! '.$metadata_path,
        );
    }

    return @errors;
}

sub execute {
    my $self = shift;
    $self->status_message('Load TCGA metadata...');

    my $metadata_path = $self->metadata_path;
    $self->status_message('Metadata path: '.$metadata_path);
    my $metadata = XMLin($self->metadata_path);
    if ( not $metadata ) {
        $self->error_message('Failed to open metadata path as XML!');
        return;
    }

    my $instrument_data_attributes = $self->_get_instrument_data_attributes_from_metadata($metadata);
    return if not $instrument_data_attributes;

    my $target_region = $self->_determine_target_region_from_library_strategy( $instrument_data_attributes->{library_strategy} );
    return if not $target_region;
    $instrument_data_attributes->{target_region} = $target_region;

    my $md5 = $self->_get_md5_from_metadata($metadata);
    $instrument_data_attributes->{bam_md5} = $md5 if $md5;

    $self->instrument_data_attributes($instrument_data_attributes);

    $self->status_message('Load TCGA metadata...done');
    return 1;
}

sub _get_instrument_data_attributes_from_metadata {
    my ($self, $metadata) = @_;

    my %attrs_and_metadata_attrs = (
        tcga_name => 'legacy_sample_id',
        import_source_name => 'center_name',
        analysis_id => 'analysis_id',
        aliquot_id => 'aliquot_id',
        participant_id => 'participant_id',
        sample_id => 'sample_id',
        library_strategy => 'library_strategy',
    );
    my %instrument_data_attributes;
    for my $attr ( keys %attrs_and_metadata_attrs ) {
        my $value = $metadata->{Result}->{ $attrs_and_metadata_attrs{$attr} };
        next if not defined $value or $value eq '';
        $instrument_data_attributes{ $attr } = $value;
    }

    return \%instrument_data_attributes;
}

sub _determine_target_region_from_library_strategy {
    my ($self, $library_strategy) = @_;

    if ( not $library_strategy ) {
        $self->status_message('No library strategy in metadata!');
        return 'none';
    }

    my $target_region;
    if ( $library_strategy eq 'WGS' or $library_strategy eq 'RNA-Seq' ) {
        $target_region = 'none';
    } elsif($library_strategy eq 'WXS') { #exome
        $target_region = 'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37';
    }
    else {
        $self->status_message('Unknown library strategy: '.$library_strategy);
        $target_region = 'none';
    }

    $self->status_message('Library strategy: '.$library_strategy);
    $self->status_message('Target region: '.$target_region);

    return $target_region;
}

sub _get_md5_from_metadata {
    my ($self, $metadata) = @_;

    Carp::confess('No metadata to retrieve MD5!') if not $metadata;

    if ( my $md5 = eval{ $metadata->{Result}->{files}->{file}->{checksum}->{content}; } ) {
        $self->status_message('MD5: '.$md5);
        return $md5;
    }

    my ($fileinfo) = eval{ grep {$_->{filename} =~ m/\.bam$/} @{$metadata->{Result}->{files}->{file}} };
    if( $fileinfo and  $fileinfo->{checksum}->{content} ) {
        $self->status_message('MD5: '.$fileinfo->{checksum}->{content});
        return $fileinfo->{checksum}->{content};
    }

    $self->warning_message("Failed to get files md5 info from metadata file");
    return;
}

1;

