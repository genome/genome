package Genome::Model::Tools::CgHub::Metadata;

use strict;
use warnings;

use Genome;

use Params::Validate ':types';
use XML::Simple;

class Genome::Model::Tools::CgHub::Metadata {
    is => 'UR::Object',
    has => {
        metadata => { 
            is => 'Hash',
            doc => 'The hash version of the metadata XML.',
        },
    },
};

sub create_from_xml {
    my ($class, $xml) = Params::Validate::validate_pos(@_, {type => SCALAR}, {type => SCALAR, called => 'XML',});

    my $metadata = XMLin($xml);
    die $class->error_message('Failed to load metadata XML!') if not $metadata;

    return $class->SUPER::create(metadata => $metadata);
}

sub create_from_file {
    my ($class, $metadata_file) = Params::Validate::validate_pos(@_, {type => SCALAR}, {type => SCALAR});

    my $xml = Genome::Sys->read_file($metadata_file);
    if ( not $xml ) {
        die $class->error_message('Failed to read in CG Query metadata XML from %s', $metadata_file);
    }

    return $class->create_from_xml($xml);
}

my %attribute_mapping = (
    uuid => 'analysis_id',
    tcga_name => 'legacy_sample_id',
    import_source_name => 'center_name',
);
sub get_attribute_value {
    my ($self, $name) = @_;

    die $self->error_message('No metadata set to get attribute value!') if not $self->metadata;
    die $self->error_message('No name to get attribute value!') if not $name;

    $name = $attribute_mapping{$name} if exists $attribute_mapping{$name};
    my $value = $self->metadata->{Result}->{$name};

    return if not defined $value or $value eq '';
    return $value;
}

sub uuid {
    return $_[0]->get_attribute_value('uuid');
}

sub _files {
    my $self = shift;

    die $self->error_message('No metadata set to get files!') if not $self->metadata;

    my $files_ref = ref $self->metadata->{Result}->{files}->{file};
    return (
        $files_ref eq 'ARRAY'
        ? @{$self->metadata->{Result}->{files}->{file}}
        : $self->metadata->{Result}->{files}->{file}
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

sub file_size_for_file_name {
    my ($self, $bam_file) = @_;
    return $self->_get_attribute_value_for_file_name($bam_file, 'filesize');
}

1;

