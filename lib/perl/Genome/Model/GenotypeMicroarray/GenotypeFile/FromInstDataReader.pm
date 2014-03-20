package Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataReader;

use strict;
use warnings;

use Genome;

sub create {
    my ($class, %params) = @_;

    my $self = bless(\%params, $class);

    if ( not $self->{instrument_data} ) {
        $class->error_message('No instrument data given to open instrument data reader!');
        return;
    }

    my $open_genotype_fh_ok = $self->_open_genotype_fh;
    return if not $open_genotype_fh_ok;

    my $resolve_headers = $self->_resolve_headers;
    return if not $resolve_headers;

    # Set the snp id mapping, if given a variation list build
    if ( $self->{variation_list_build} ) {
        $self->{snp_id_mapping} = Genome::InstrumentData::Microarray->get_snpid_hash_for_variant_list(
            $self->{instrument_data}, $self->{variation_list_build}
        );
    }
    $self->{snp_id_mapping} ||= {};

    return $self;
}

sub _open_genotype_fh {
    my $self = shift;

    my $genotype_file = eval{ $self->{instrument_data}->attributes(attribute_label => 'genotype_file')->attribute_value; };
    if ( not $genotype_file or not -s $genotype_file ) {
        $self->error_message('No genotype file for instrument data! '.Data::Dumper::Dumper($self->{instrument_data}));
        return;
    }

    my $fh = eval { Genome::Sys->open_file_for_reading($genotype_file) };
    if (!$fh or $@) {
        $self->error_message("Can't open file $genotype_file for reading: $@");
        return;
    }

    $self->{genotype_file} = $genotype_file;
    $self->{genotype_fh} = $fh;

    return 1;
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
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split(',', $header_line);
    $self->{headers} = \@headers;

    return 1;
}

sub read {
    my $self = shift;

    my $line = $self->{genotype_fh}->getline;
    return if not $line;
    chomp $line;

    my %genotype;
    @genotype{ @{$self->{headers}} } = split(',', $line);

    # The id is from the snp mapping or the genotype's snp_name
    $genotype{id} = delete $genotype{snp_name};
    if ( exists $self->{snp_id_mapping}->{ $genotype{id} }) {
        $genotype{id} = $self->{snp_id_mapping}->{ $genotype{id} };
    }

    return \%genotype;
}

1;

