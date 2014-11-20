package Genome::InstrumentData::Command::Import::WorkFlow::Tcga::Metadata;

use strict;
use warnings;

use Genome;

require File::Spec;
require File::Temp;
use XML::Simple;

class Genome::InstrumentData::Command::Import::WorkFlow::Tcga::Metadata { 
    is => 'UR::Object',
    has => {
        metadata_file => {
            is => 'Text',
            doc => 'The full path of a XML metadata file containing information about a particular TCGA run including sequence file information. A temp file will be used if not given.',
        },
    },
    has_optional => {
        uuid => {
            is => 'Text',
            doc => 'The UUID (analysis_id) for the TCGA run used to retrieve the metadata file.',
        },
    },
    has_optional_transient => {
        _metadata => { 
            is => 'Hash',
            doc => 'The hash version of the metadata XML.',
        },
    },
};

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(%params);
    return if not $self;

    my $retrieve = $self->_retrieve_metadata_file;
    return if not $retrieve;

    my $load = $self->_load_metadata_file;
    return if not $load;

    return $self;
}

sub _retrieve_metadata_file {
    my $self = shift;

    my $metadata_file = $self->metadata_file;
    return 1 if defined $metadata_file and -s $self->metadata_file;

    my $uuid = $self->uuid;
    die $self->error_message('Need UUID or existing metadata file!') if not defined $uuid;

    my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    $metadata_file = $self->metadata_file( File::Spec->catfile($tmp_dir, 'metadata.xml') ) if not defined $metadata_file;
    my $output_file = File::Spec->catfile($tmp_dir, 'cgquery.out');

    local $ENV{PATH} = $ENV{PATH} . ':/cghub/bin';
    my $cmd = "cgquery -o $metadata_file analysis_id=$uuid > $output_file";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $output_file ) {
        $self->error_message($@) if $@;
        die $self->error_message('Failed to run cgquery to retrieve metadata file!');
    }

    my $output_fh = Genome::Sys->open_file_for_reading($output_file);
    my $number_of_hits = 0;
    while ( my $line = $output_fh->getline ) {
        if ( $line =~ /Matching Objects\s+:\s+(\d+)/ ) {
            $number_of_hits = $1;
            last;
        }
    }
    $output_fh->close;
    if ( $number_of_hits == 0 ) {
        die $self->error_message("Failed to find uuid ($uuid) on CG Hub!");
    }


    return 1;
}

sub _load_metadata_file {
    my $self = shift;

    my $metadata_file = $self->metadata_file;
    Genome::Sys->validate_file_for_reading($metadata_file);
    my $metadata = XMLin($metadata_file);
    if ( not $metadata ) {
        die $self->error_message('Failed to open metadata file! '.$metadata_file);
    }

    $self->_metadata($metadata);
    $self->uuid( $metadata->{Result}->{analysis_id} ) if not defined $self->uuid;

    return 1;
}

my %attribute_mapping = (
    uuid => 'analysis_id',
    tcga_name => 'legacy_sample_id',
    import_source_name => 'center_name',
);
sub get_attribute_value {
    my ($self, $name) = @_;

    die $self->error_message('No metadata set to get attribute value!') if not $self->_metadata;
    die $self->error_message('No name to get attribute value!') if not $name;

    if ( $self->can($name) ) {
        return $self->$name;
    }

    $name = $attribute_mapping{$name} if exists $attribute_mapping{$name};
    my $value = $self->_metadata->{Result}->{$name};

    return if not defined $value or $value eq '';
    return $value;
}

sub reference_assembly_shortname {
    my $self = shift;

    my $reference_shortname = $self->get_attribute_value('refassem_short_name');
    return $reference_shortname if $reference_shortname;

    $reference_shortname = $self->_metadata->{Result}->{analysis_xml}->{ANALYSIS_SET}->{ANALYSIS}->{ANALYSIS_TYPE}->{REFERENCE_ALIGNMENT}->{ASSEMBLY}->{STANDARD}->{shortname};
    return $reference_shortname if $reference_shortname;

    return;
}

sub reference_assembly_version {
    my $self = shift;

    my $reference_assembly_shortname = $self->reference_assembly_shortname;
    return if not $reference_assembly_shortname;

    if ( $reference_assembly_shortname =~ /hg(\d\d)/i ) {
        return $1 + 18;
    }
    elsif ( $reference_assembly_shortname =~ /ncbi(\d\d)/i ) {
        return $1;
    }

    return;
}

sub target_region {
    # This will need to be updated
    my $self = shift;

    die $self->error_message('No metadata set to get target region!') if not $self->_metadata;

    my $target_region;
    my $library_strategy = $self->get_attribute_value('library_strategy');
    die $self->error_message('No library strategy in metadata to resolve target region!') if not $library_strategy;
    if ( $library_strategy eq 'WGS' or $library_strategy eq 'RNA-Seq' ) {
        $target_region = 'none';
    } elsif($library_strategy eq 'WXS') { #exome
        my $reference_assembly_version = $self->reference_assembly_version;
        $target_region = ( $reference_assembly_version eq '37' )
        ? 'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37'
        : 'agilent sureselect exome version 2 broad refseq cds only';
    }
    else {
        # unknown library strategy
        $target_region = 'none';
    }

    return $target_region;
}

sub _files {
    my $self = shift;

    die $self->error_message('No metadata set to get files!') if not $self->_metadata;

    my $files_ref = ref $self->_metadata->{Result}->{files}->{file};
    return (
        $files_ref eq 'ARRAY'
        ? @{$self->_metadata->{Result}->{files}->{file}}
        : $self->_metadata->{Result}->{files}->{file}
    );
}

sub file_names {
    my $self = shift;
    my @files = $self->_files;
    return map { $_->{filename} } @files;
}
sub bam_file_names {
    my $self = shift;
    my @file_names = $self->file_names;
    my @bam_file_names = grep { $_ =~ m/\.bam$/} @file_names;
    return wantarray ? @bam_file_names : $bam_file_names[0];
}

sub _get_attribute_value_for_file_name {
    my ($self, $file_name, @attr_names) = @_;

    die $self->error_message('No file name given to get attribute value!') if not $file_name;
    die $self->error_message('No attribute name(s) to get attribute value!') if not @attr_names;

    my @files = $self->_files;
    return if not @files;

    my @matching_files = grep { $_->{filename} eq $file_name } @files;
    return if not @matching_files;

    my @values;
    if ( @attr_names == 1 ) {
        @values = map { $_->{$attr_names[0]} } @matching_files;
    }
    elsif ( @attr_names == 2 ) { 
        @values = map { $_->{$attr_names[0]}->{$attr_names[1]} } @matching_files;
    }
    else {
        die $self->error_message('Unknown number of attribute names given to _get_attribute_value_for_file_name');
    }

    return if not @values;
    die $self->error_message("Got multiple values (@values) for $file_name @attr_names!") if @values == 2;

    return $values[0];
}

sub checksum_type_for_file_name {
    my ($self, $file_name) = @_;
    my $checksum_type = $self->_get_attribute_value_for_file_name($file_name, 'checksum', 'type');
    return uc $checksum_type;
}

sub checksum_content_for_file_name {
    my ($self, $file_name) = @_;
    return $self->_get_attribute_value_for_file_name($file_name, 'checksum', 'content');
}

sub filesize_in_kb_for_file_name {
    my ($self, $bam_file) = @_;
    my $filesize = $self->_get_attribute_value_for_file_name($bam_file, 'filesize');
    return if not defined $filesize;
    return sprintf('%d', ($filesize / 1024));
}

1;

