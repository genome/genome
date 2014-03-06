package Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate { 
    is => 'UR::Object',
    has => {
        input => { is => 'Text', },
        variation_list_build => { is => 'Genome::Model::Build::ImportedVariationList', },
        snp_id_mapping => { is => 'Hash', },
        _genotype_fh => { is => 'IO::File', },
        _headers => { is => 'Array', },
        _sample => { is => 'Genome::Sample', },
        _genotypes => { is => 'Hash', default_value => {}, },
        _order => { is => 'Array', },
        _position => { is => 'Integer', default_value => 0, },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $open_ok = $self->_open_genotype_file;
    return if not $open_ok;

    my $headers_ok = $self->_resolve_headers;
    return if not $headers_ok;

    my $sample_ok = $self->_resolve_sample;
    return if not $sample_ok;

    my $load_genotypes = $self->_load_genotypes;
    return if not $load_genotypes;

    my $annotate_genotypes = $self->_annotate_genotypes;
    return if not $annotate_genotypes;

    return $self;
}

BEGIN {
    *next = \&read;
}
sub read {
    my $self = shift;

    my $position = $self->_position;
    my $id;
    do {
        $id = $self->_order->[$position++];
    } until not $id or not $self->_genotypes->{$id}->{ignored};
    $self->_position($position);

    return if not defined $id;
    return $self->_genotypes->{$id};
}

sub _open_genotype_file {
    my $self = shift;

    my $genotype_file = $self->input;
    if ( not -s $genotype_file ) {
        $self->error_message('Genotype file file does not exist! '.$genotype_file);
        return;
    }

    my $genotype_fh = eval{ Genome::Sys->open_file_for_reading($genotype_file); };
    if ( not $genotype_fh ) {
        $self->error_message("Failed to open reader for genotype file: $genotype_file): $@");
        return;
    }

    $self->_genotype_fh($genotype_fh);

    return 1;
}

sub _resolve_headers {
    my $self = shift;

    my $header_line;
    my $genotype_fh = $self->_genotype_fh;
    do { $header_line = $genotype_fh->getline; } until not $header_line or $header_line =~ /,/;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }

    chomp $header_line;
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split(',', $header_line);
    $self->_headers(\@headers);

    return 1;
}

sub _resolve_sample {
    my $self = shift;

    my $genotype = $self->_load_genotype;
    if ( not $genotype ) {
        $self->error_message('No genotype found after header!');
        return;
    }

    my %sample_params;
    my $sample_id = $genotype->{sample_id};
    if ( defined $genotype->{sample_id} ) {
        $sample_params{id} = $genotype->{sample_id};
    }
    elsif ( defined $genotype->{sample_name} ) {
        $sample_params{name} = $genotype->{sample_name};
    }
    else {
        $self->error_message('No sample id or name in genotype!');
        return;
    }

    my $sample = Genome::Sample->get(%sample_params);
    if ( not $sample ) {
        $self->error_message('No sample for params! '.Data::Dumper::Dumper(\%sample_params));
        return;
    }
    $self->_sample($sample);

    return 1;
}

sub _load_genotype {
    my $self = shift;

    my $line = $self->_genotype_fh->getline;
    return if not $line;

    chomp $line;
    my %genotype;
    @genotype{@{$self->_headers}} = split(',', $line);

    # The id is from the snp mapping or the genotype's snp_name
    if ( $self->snp_id_mapping and exists $self->snp_id_mapping->{ $genotype{snp_name} }) {
        $genotype{id} = $self->snp_id_mapping->{ delete $genotype{snp_name} };
    } else {
        $genotype{id} = delete $genotype{snp_name};
        $genotype{id} =~ s/^(rs\d+)\D*$/$1/; #borrowed from GSC::Genotyping::normalize_to
    }

    if ( exists $self->_genotypes->{ $genotype{id} } ) {
        Carp::confess( $self->error_message('Already have a genotype for snp id: '.Dumper(\%genotype, $self->genotypes->{ $genotype{id} })) );
    }

    delete $genotype{'chr'};
    $genotype{alleles} = $genotype{allele1}.$genotype{allele2};

    $self->_genotypes->{ $genotype{id} } = \%genotype;

    return \%genotype;
}

sub _load_genotypes {
    my $self = shift;

    my $snp_id_mapping = $self->snp_id_mapping;

    my $genotype_fh = $self->_genotype_fh;
    my @headers = @{$self->_headers};
    my $genotypes = $self->_genotypes;

    my $genotype;
    do {
        $genotype = $self->_load_genotype;
    } while $genotype;

    if ( not %{$self->_genotypes} ) {
        $self->error_message("No genotypes found in genotype file! ".$self->get_original_input);
        return;
    }

    return 1;
}

sub _annotate_genotypes {
    my $self = shift;

    my $genotypes = $self->_genotypes;
    Carp::confess('No genotypes!') if not $genotypes or not %$genotypes;

    my $variation_list_build = $self->variation_list_build;
    my $snvs_file = $variation_list_build->snvs_bed;
    if ( not $snvs_file ) {
        $self->error_message('No snvs file (snvs_bed) for build: '.$variation_list_build->__display_name__);
        return;
    }

    my $dbsnp_fh = eval{ Genome::Sys->open_file_for_reading($snvs_file); };
    if ( not $dbsnp_fh ) {
        $self->error_message("Failed to open file: $snvs_file");
        return;
    }

    my $variant_id_pos = ( $variation_list_build->version and $variation_list_build->version eq 130 ? 8 : 7 );
    my $reference_sequence_build = $variation_list_build->reference;
    my %order;
    my $cnt = 0;
    while ( my $line = $dbsnp_fh->getline ) {
        chomp $line;
        my @tokens = split(/\s+/, $line);
        my $variant_id = $tokens[$variant_id_pos];
        my $genotype = $genotypes->{$variant_id};
        next if not $genotype;

        if ( exists $order{$variant_id} ) {
            if ( $genotype->{position} != $tokens[0] or $genotype->{position} != $tokens[2] ) {
                $genotype->{ignored} = 1;
            }
            next;
        }

        $genotype->{chromosome} = $tokens[0];
        $genotype->{position} = $tokens[2];
        $genotype->{reference} = $reference_sequence_build->sequence(
            $genotype->{chromosome}, $genotype->{position}, $genotype->{position}
        );

        # Order
        $order{$variant_id} = $cnt++;
    }

    my @order = sort { $order{$a} <=> $order{$b} } grep { !$genotypes->{$_}->{ignored} } keys %order;
    if ( not @order ) {
        $self->error_message("All genotypes are duplicates in variant list! ".$snvs_file);
        return;
    }
    $self->_order(\@order);

    return 1;
}

1;

