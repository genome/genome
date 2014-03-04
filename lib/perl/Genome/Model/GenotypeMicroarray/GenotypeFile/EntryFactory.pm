package Genome::Model::GenotypeMicroarray::GenotypeFile::EntryFactory;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Entry;
use Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper;

class Genome::Model::GenotypeMicroarray::GenotypeFile::EntryFactory {
    has => {
        sample_name => { is => 'Text', },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my @errors = $self->__errors__;
    if ( @errors ) {
        $self->error_message( $errors[0]->desc );
        return;
    }

    return $self;
}

sub build_entry {
    my ($self, $genotype) = @_;

    Carp::confess('No genotype!') if not $genotype;

    my $type = ref($genotype);
    if ( not $type ) {
        Carp::confess('Need genotype as HASH or VCF enrtry! ');
    }

    #TODO verify sample_name

    if ( $type eq 'Genome::File::Vcf::Entry' ) {
        return $self->_build_from_vcf_entry($genotype);
    }
    elsif ( $type eq 'HASH' ) {
        return $self->_build_from_hash($genotype);
    }
    else {
        Carp::confess('Unknown genotype type! '.$type);
    }

}

#< From VCF >#
sub _build_from_vcf_entry {
    my ($self, $genotype) = @_;

    $genotype->{sample_name} = $self->sample_name;
    $genotype->{id} = $genotype->{identifiers}->[0];
    $genotype->{chromosome} = $genotype->{chrom};
    $genotype->{reference} = $genotype->{reference_allele};
    if ( not $genotype->{_sample_data} ) {
        $genotype->sample_data;
    }
    $genotype->{alleles} = $genotype->{_sample_data}->[0]->[1];
    my @alleles = split(//, $genotype->{alleles});
    $genotype->{allele1} = $alleles[0];
    $genotype->{allele2} = $alleles[1];

    for my $info_field ( @{Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->supported_info_fields} ) {
        my $value =  $genotype->info( $info_field->{name});
        $genotype->{$info_field->{name}} = defined $value ? $value : 'NA';
    }

    return $genotype;
}

#< From HASH >#
sub _build_from_hash {
    my ($self, $genotype) = @_;

    $genotype->{sample_name} = $self->sample_name;
    $genotype->{identifiers} = [ $genotype->{id} ];
    $genotype->{chrom} = $genotype->{chromosome};
    $genotype->{reference_allele} = $genotype->{reference};
    $genotype->{alternate_alleles} = $self->alts_for_genotype($genotype);
    $genotype->{quality} = '.';
    $genotype->{_filter} = [];
    $genotype->{info_fields} = $self->info_hash_for_genotype($genotype),
    $genotype->{_format} = Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->formats;
    $genotype->{_format_key_to_idx} = Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->format_key_idx;
    $genotype->{_sample_data} = $self->sample_data_for_genotype($genotype);

    return bless($genotype, 'Genome::File::Vcf::Entry');
}

sub alts_for_genotype {
    my ($self, $genotype) = @_;

    # alleles are dashes
    return '.' if $genotype->{allele1} eq '-' and $genotype->{allele2} eq '-';

    # see which match the ref
    my @alts;
    for my $allele_key ( map{ 'allele'.$_ } (1..2) ) {
        next if $genotype->{$allele_key} eq '-';
        push @alts, $genotype->{$allele_key} if $genotype->{$allele_key} ne $genotype->{reference_allele};
    }

    # both match!
    return [] if not @alts;
    # retrun alts
    return [ sort { $a cmp $b } @alts ];
}

sub sample_data_for_genotype {
    my ($self, $genotype) = @_;

    my $allele1_matches_ref = $genotype->{allele1} eq $genotype->{reference};
    my $allele2_matches_ref = $genotype->{allele2} eq $genotype->{reference};

    my $gt;
    if ( $genotype->{allele1} eq '-' and $genotype->{allele2} eq '-' ) {
        $gt = './.';
    }
    elsif ( $allele1_matches_ref and $allele2_matches_ref ) { # ref and alleles match
        $gt = '0/0'; # homozygous ref
    }
    elsif ( $allele1_matches_ref or $allele2_matches_ref ) { # ref matches one allele, but alleles are not the same
        $gt = '0/1'; # heterozygous ref
    }
    elsif ( $genotype->{allele1} eq $genotype->{allele2} ) { # alleles match, but not the ref
        $gt = '1/1'; # homozygous alt
    }
    else { # not $allele1_matches_ref and not $allele2_matches_ref
        $gt = '1/2';
    }

    return [ [ $gt, $genotype->{alleles} ] ];

}

sub info_hash_for_genotype {
    my ($self, $genotype) = @_;

    my (%info_hash, @order);
    for my $info_field ( @{Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->supported_info_fields} ) { 
        my $value = $genotype->{ $info_field->{name} };
        next if not defined $value or $value eq 'NA'; # FIXME convert undef/NA to . ??
        $info_hash{ $info_field->{id} } = $value;
        push @order, $info_field->{id};
    }

    return { hash => \%info_hash, order => \@order, };
}

1;

