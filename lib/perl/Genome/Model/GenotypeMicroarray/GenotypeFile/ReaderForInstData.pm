package Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForInstData;

use strict;
use warnings;

use parent 'Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForBase';

sub source_type {
    return 'instrument_data';
}

sub get_genotype_file {
    my $self = shift;
    return $self->instrument_data->genotype_file;
}

sub variation_list_build {
    return $_[0]->{variation_list_build};
}

sub instrument_data {
    return $_[0]->{instrument_data};
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    if ( $self->variation_list_build ) {
        $self->{snp_id_mapping} = Genome::InstrumentData::Microarray->get_snpid_hash_for_variant_list(
            $self->instrument_data, $self->variation_list_build
        );
    }
    $self->{snp_id_mapping} ||= {};

    return $self;
}

sub _resolve_headers {
    my $self = shift;

    my $header_line;
    do { $header_line = $self->{genotype_fh}->getline; } until not $header_line or $header_line =~ /,/;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }

    chomp $header_line;
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split($self->{delimiter}, $header_line);
    $self->{headers} = \@headers;

    return 1;
}

sub read {
    my $self = shift;

    my $genotype = $self->SUPER::read;
    return if not $genotype;

    $genotype->{id} = delete $genotype->{snp_name}; # set snp name as id
    if ( exists $self->{snp_id_mapping}->{ $genotype->{id} }) { # get corrected id
        $genotype->{id} = $self->{snp_id_mapping}->{ $genotype->{id} };
    }

    return $genotype;
}

1;

