package Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataReader;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataReader {
    is => 'Genome::Utility::IO::SeparatedValueReader',
    has_optional => {
        variation_list_build => { is => 'Genome::Model::Build::ImportedVariationList', },
        snp_id_mapping => { is => 'Hash', },
    },
};

sub create {
    my ($class, %params) = @_;

    my $instrument_data = delete $params{instrument_data};
    if ( not $instrument_data ) {
        $class->error_message('No instrument data given to open instruemnt data reader!');
        return;
    }

    my $genotype_file = eval{ $instrument_data->attributes(attribute_label => 'genotype_file')->attribute_value; };
    if ( not $genotype_file or not -s $genotype_file ) {
        $class->error_message('No genotype file for instrument data! '.Data::Dumper::Dumper($instrument_data));
        return;
    }

    $params{input} = $genotype_file;

    $params{headers} = []; # holder

    my $self = $class->SUPER::create(%params);
    return if not $self;

    my $resolve_headers = $self->_resolve_headers;
    return if not $resolve_headers;

    # Set the snp id mapping, if given a variation list build
    if ( $self->variation_list_build ) {
        $self->snp_id_mapping( 
            Genome::InstrumentData::Microarray->get_snpid_hash_for_variant_list(
                $instrument_data, $self->variation_list_build
            )
        );
    }


    return $self;
}

sub _resolve_headers {
    my $self = shift;

    my $header_line;
    do { $header_line = $self->getline; } until not $header_line or $header_line =~ /,/;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }

    chomp $header_line;
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split(',', $header_line);
    $self->headers(\@headers);

    return 1;
}

sub read {
    my $self = shift;

    # Get next genotype
    my $genotype = $self->next;
    return if not $genotype;

    # The id is from the snp mapping or the genotype's snp_name
    if ( $self->snp_id_mapping and exists $self->snp_id_mapping->{ $genotype->{snp_name} }) {
        $genotype->{id} = $self->snp_id_mapping->{ delete $genotype->{snp_name} };
    } else {
        $genotype->{id} = delete $genotype->{snp_name};
        $genotype->{id} =~ s/^(rs\d+)\D*$/$1/; #borrowed from GSC::Genotyping::normalize_to
    }

    return $genotype;
}

1;

