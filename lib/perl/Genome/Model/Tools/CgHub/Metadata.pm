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
            default_value => {},
            doc => 'The hash version of the metadata XML.',
        },
    },
};

sub add_file {
    my ($self, $metadata_file) = Params::Validate::validate_pos(@_, {type => OBJECT}, {type => SCALAR, called => 'file'});

    my $xml = Genome::Sys->read_file($metadata_file);
    if ( not $xml ) {
        die $self->error_message('Failed to read in CG Query metadata XML from %s', $metadata_file);
    }

    return $self->add_xml($xml);
}

sub add_xml {
    my ($self, $xml) = Params::Validate::validate_pos(@_, {type => OBJECT}, {type => SCALAR, called => 'XML',});

    my $metadata = XMLin($xml);
    die $self->error_message('Failed to load metadata XML!') if not $metadata;
    $self->_add_metadata($metadata);

    return 1;
}

sub _add_metadata {
    my ($self, $incoming_metadata) = Params::Validate::validate_pos(@_, {type => OBJECT}, {type => HASHREF, called => 'metadata',});

    my $results;
    if ( not $incoming_metadata->{Hits} or $incoming_metadata->{Hits} == 0) {
        return 1;
    }
    elsif ( $incoming_metadata->{Hits} == 1) {
        $results = [ $incoming_metadata->{Result} ];
    }
    else {
        $results = [ values %{$incoming_metadata->{Result}} ];
    }

    my $metadata = $self->metadata;
    for my $result ( @$results ) {
        for my $id_by (qw/ analysis_id legacy_sample_id /) {
            $metadata->{ $result->{$id_by} } = $result;
        }
    }
    $self->metadata($metadata);

    return 1;
}

sub get_attribute_value {
    my ($self, $lookup_id, $name) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => SCALAR},
    );

    my $result = $self->metadata->{$lookup_id};
    return if not $result;
    my $value = $result->{$name};

    return if not defined $value or $value eq '';
    return $value;
}

sub _files {
    my ($self, $lookup_id) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => SCALAR});

    my $result = $self->metadata->{$lookup_id};
    return if not $result;

    my $files = $result->{files}->{file};
    my $files_ref = ref $files;
    return ( $files_ref eq 'ARRAY' ? @$files : $files );
}

sub file_names {
    my $self = shift;
    my @files = $self->_files(@_);
    return map { $_->{filename} } @files;
}

sub bam_file_names {
    my $self = shift;
    my @file_names = $self->file_names(@_);
    my @bam_file_names = grep { $_ =~ m/\.bam$/} @file_names;
    return wantarray ? @bam_file_names : $bam_file_names[0];
}

sub _get_attribute_value_for_file_name {
    my ($self, $lookup_id, $file_name, $attr_names) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => SCALAR}, {type => ARRAYREF},
    );

    my @files = $self->_files($lookup_id);
    return if not @files;

    my @matching_files = grep { $_->{filename} eq $file_name } @files;
    return if not @matching_files;

    my @values;
    if ( @$attr_names == 1 ) {
        @values = map { $_->{$attr_names->[0]} } @matching_files;
    }
    elsif ( @$attr_names == 2 ) { 
        @values = map { $_->{$attr_names->[0]}->{$attr_names->[1]} } @matching_files;
    }
    else {
        die $self->error_message('Unknown number of attribute names given to _get_attribute_value_for_file_name');
    }

    return if not @values;
    die $self->error_message("Got multiple values (%s) for %s %s", join(' ', @values), $file_name, join(' ', @$attr_names)) if @values == 2;

    return $values[0];
}

sub checksum_type_for_file_name {
    my ($self, $lookup_id, $file_name) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => SCALAR},
    );
    my $checksum_type = $self->_get_attribute_value_for_file_name($lookup_id, $file_name, [qw/ checksum type /]);
    return uc $checksum_type;
}

sub checksum_content_for_file_name {
    my ($self, $lookup_id, $file_name) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => SCALAR},
    );
    return $self->_get_attribute_value_for_file_name($lookup_id, $file_name, [qw/ checksum content /]);
}

sub file_size_for_file_name {
    my ($self, $lookup_id, $bam_file) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => SCALAR},
    );
    return $self->_get_attribute_value_for_file_name($lookup_id, $bam_file, [qw/ filesize /]);
}

1;

