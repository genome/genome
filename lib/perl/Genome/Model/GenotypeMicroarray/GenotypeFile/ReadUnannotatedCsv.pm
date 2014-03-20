package Genome::Model::GenotypeMicroarray::GenotypeFile::ReadUnannotatedCsv;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Reader;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReadUnannotatedCsv { 
    is => 'Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataReader',
    has => {
        _vcf_reader => { is => 'Genome::File::Vcf::Reader', },
        _genotypes => { is => 'Hash', default_value => {}, },
        _order => { is => 'Array', },
        _position => { is => 'Integer', default_value => 0, },
    },
};

sub header {
    my ($self, $header) = @_;

    if ( $header ) {
        $self->_vcf_reader->header($header);
    }

    return $self->_vcf_reader->header;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $open_vcf_reader = $self->_open_vcf_reader;
    return if not $open_vcf_reader;

    return $self;
}

BEGIN {
    *next = \&read;
}

sub read {
    my $self = shift;

    if ( not defined $self->_order ) {
        my $load_genotypes = $self->_load_genotypes;
        die if not $load_genotypes;

        my $annotate_genotypes = $self->_annotate_genotypes;
        die if not $annotate_genotypes;
    }

    while ( my $variant_id = shift @{$self->_order} ) {
        my $entry = $self->_create_vcf_entry( $self->_genotypes->{$variant_id} );
        return $entry if $entry;
    }

    return;
}

sub _create_vcf_entry {
    my ($self, $genotype) = @_;

    my $entry = Genome::File::Vcf::Entry->new($self->header, $genotype->{line});

    # Skip INDELs
    return if $entry->has_indel;

    # Add GT to genotype
    $genotype->{genotype} = $self->_gt_for_genotype($genotype, $entry);

    # Add genotype data to entry
    for my $field ( Genome::Model::GenotypeMicroarray->format_types ) {
        $entry->add_format_field($field->{id});
        my $value = $genotype->{ $field->{name} };
        $value = '.' if not defined $value or $value eq '';
        $entry->set_sample_field(0, $field->{id}, $value);
    }

    return $entry;
}

sub _open_vcf_reader {
    my $self = shift;

    my $variation_list_build = $self->variation_list_build;
    my $snvs_vcf = $variation_list_build->snvs_vcf;
    if ( not $snvs_vcf or not -s $snvs_vcf ) {
        $self->error_message('No SNVs VCF for variation list  build! '.$variation_list_build->__display_name__);
        return;
    }

    my $vcf_reader = eval{ Genome::File::Vcf::Reader->new($snvs_vcf); };
    if ( not $vcf_reader ) {
        $self->error_message("Failed to open SNVs VCF file! $snvs_vcf");
        return;
    }
    $self->_vcf_reader($vcf_reader);

    return 1;
}

sub _load_genotype {
    my $self = shift;

    my $genotype = $self->SUPER::read;
    return if not $genotype;

    if ( exists $self->_genotypes->{ $genotype->{id} } ) {
        Carp::confess( $self->error_message('Already have a genotype for snp id: '.Data::Dumper::Dumper($genotype, $self->genotypes->{ $genotype->{id} })) );
    }

    delete $genotype->{'chr'};
    $genotype->{alleles} = $genotype->{allele1}.$genotype->{allele2};

    $self->_genotypes->{ $genotype->{id} } = $genotype;

    return $genotype;
}

sub _load_genotypes {
    my $self = shift;

    my $snp_id_mapping = $self->snp_id_mapping;

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

    my $vcf_reader = $self->_vcf_reader;
    my $cnt = 0;
    while ( my $line = $vcf_reader->_getline ) {
        my ($variant_id) = split(',', (split(/\t/, $line))[2]);
        my $genotype = $genotypes->{$variant_id};

        # Skip if not in variation list
        next if not $genotype;

        if ( $genotypes->{$variant_id}->{seen} ) {
            # remove order value if seen
            delete $genotypes->{$variant_id}->{order};
        }
        else {
            # set order
            $genotypes->{$variant_id}->{order} = ++$cnt;
        }
        # increment seen
        $genotypes->{$variant_id}->{seen}++;

        # Save line to create vcf entry
        chomp $line;
        $genotypes->{$variant_id}->{line} = $line;
    }

    my @order = map { $_->{id} } sort { $a->{order} <=> $b->{order} } grep { exists $_->{order} } values %$genotypes;
    $self->_order(\@order);
    if ( not @order ) {
        $self->error_message("All genotypes are duplicates in variant list! ".$vcf_reader->{name});
        return;
    }

    return 1;
}

sub _gt_for_genotype {
    my ($self, $genotype, $entry) = @_;

    my %alleles_idx = (
        '-' => '.',
        $entry->{reference_allele} => 0,
    );
    @alleles_idx{ @{$entry->{alternate_alleles}} } = ( 1..@{$entry->{alternate_alleles}} );

    my @gt_idx;
    for my $allele ( map { $genotype->{$_} } (qw/ allele1 allele2 /) ) { 
        if ( not exists $alleles_idx{$allele} ) {
            push @{$entry->{alternate_alleles}}, $allele;
            $alleles_idx{$allele} = scalar(@{$entry->{alternate_alleles}});
        }
        push @gt_idx, $alleles_idx{$allele};
    }

    return join('/', @gt_idx);
}

1;

