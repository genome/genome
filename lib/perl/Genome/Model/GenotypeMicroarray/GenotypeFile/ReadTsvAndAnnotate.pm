package Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsvAndAnnotate { 
    is => 'UR::Object',
    has => {
        input => { is => 'Text', },
        separator => { is => 'Text', value => "\t", },
        variation_list_build => { is => 'Genome::Model::Build::ImportedVariationList', },
        snp_id_mapping => { is => 'Hash', },
        _genotypes => { is => 'Hash', default_value => {}, },
        _order => { is => 'Array', },
        _position => { is => 'Integer', default_value => 0, },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $load_genotypes = $self->_load_genotypes;
    return if not $load_genotypes;

    my $annotate_genotypes = $self->_annotate_genotypes;
    return if not $annotate_genotypes;

    return $self;
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

sub _load_genotypes {
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

    my $snp_id_mapping = $self->snp_id_mapping;
    my $header_line;
    do { $header_line = $genotype_fh->getline; } until not $header_line or $header_line =~ /,/;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }
    chomp $header_line;
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split(',', $header_line);

    my $genotypes = $self->_genotypes;
    while ( my $line = $genotype_fh->getline ) {
        chomp $line;
        my %genotype;
        @genotype{@headers} = split(',', $line);
        # The id is from the snp mapping or the genotype's snp_name
        if($snp_id_mapping and exists $snp_id_mapping->{ $genotype{snp_name} }) {
            $genotype{id} = $snp_id_mapping->{ $genotype{snp_name} };
        } else {
            $genotype{id} = $genotype{snp_name};
            $genotype{id} =~ s/^(rs\d+)\D*$/$1/; #borrowed from GSC::Genotyping::normalize_to
        }

        if ( exists $genotypes->{ $genotype{id} } ) {
            $self->error_message('Already have a genotype for snp id: '.Dumper(\%genotype, $genotypes->{ $genotype{id} }));
            return;
        }
        $genotype{alleles} = $genotype{allele1}.$genotype{allele2};
        $genotypes->{ $genotype{id} } = \%genotype;
    }

    if ( not %$genotypes ) {
        $self->error_message("No genotypes found in input! ".$self->get_original_input);
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
    my @order;
    my $ignored = 0;
    while ( my $line = $dbsnp_fh->getline ) {
        chomp $line;
        my @tokens = split(/\s+/, $line);
        my $variant_id = $tokens[$variant_id_pos];
        my $genotype = $genotypes->{$variant_id};
        next if not $genotype;

        if ( exists $genotype->{annotated} ) {
            if ( $genotype->{position} != $tokens[2] ) {
                $ignored++;
                $genotype->{ignored} = 1;
            }
            next;
        }

        $genotype->{chromosome} = $tokens[0];
        $genotype->{position} = $tokens[2];
        $genotype->{'ref'} = $reference_sequence_build->sequence(
            $genotype->{chromosome}, $genotype->{position}, $genotype->{position}
        );
        $genotype->{annotated} = 1;
        push @order, $variant_id;
    }

    $self->_order(\@order);
    if ( not @order ) {
        $self->error_message("All genotypes are duplicates in variant list! ".$snvs_file);
        return;
    }

    return 1;
}

1;

