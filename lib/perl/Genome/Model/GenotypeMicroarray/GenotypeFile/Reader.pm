package Genome::Model::GenotypeMicroarray::GenotypeFile::Reader;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Entry;

class Genome::Model::GenotypeMicroarray::GenotypeFile::Reader { 
    is => 'UR::Object',
    has => {
        reader => { is => 'UR::Object', },
    },
    has_optional_transient => {
        _read_from_reader => { is => 'Code', },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my @errors = $self->__errors__;
    if ( @errors ) {
        Carp::confess( join ("\n", map { $_->__display_name__ } @errors) );
    }

    $self->_resolve_read_sub;

    return $self;
}

sub read {
    return $_[0]->_read_from_reader->();
}

sub _resolve_read_sub {
    my $self = shift;

    my $reader_class = $self->reader->class;
    if ( $reader_class eq 'Genome::File::Vcf::Reader' ) {
        $self->_read_from_reader(
            sub {
                return $self->reader->next;
            }
        );
    }
    elsif ( $reader_class =~ /^Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsv/ ) {
        $self->_read_from_reader(
            sub {
                my $genotype = $self->reader->read;
                return if not $genotype;
                return $self->_convert_genotype_hash_to_vcf_entry($genotype);
            }
        );
    }
    else {
        Carp::confess('Unknown reader! '.$self->reader);
    }

}

sub _convert_genotype_hash_to_vcf_entry {
    my ($self, $genotype) = @_;

    $genotype->{identifiers} = [ $genotype->{id} ];
    $genotype->{chrom} = $genotype->{chromosome};
    $genotype->{reference_allele} = $genotype->{reference};
    $genotype->{alternate_alleles} = $self->_alts_for_genotype($genotype);
    $genotype->{quality} = '.';
    $genotype->{_filter} = [];
    $genotype->{info_fields} = {
        hash => { 
            map { $_->{id} => $genotype->{ $_->{name} } } # FIXME convert NA to . ??
            grep { defined $genotype->{$_->{name}} and $genotype->{$_->{name}} ne 'NA' }
            @{Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->supported_info_fields}, 
        },
        order =>Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->info_order,
    };
    $genotype->{_format} = [ 'GT' ];
    $genotype->{_sample_data} = [ [ $self->_format_for_genotype($genotype) ] ];

    return bless $genotype, 'Genome::File::Vcf::Entry';
}

sub _alts_for_genotype {
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

sub _format_for_genotype {
    my ($self, $genotype) = @_;

    my $allele1_matches_ref = $genotype->{allele1} eq $genotype->{reference};
    my $allele2_matches_ref = $genotype->{allele2} eq $genotype->{reference};

    if ( $genotype->{allele1} eq '-' and $genotype->{allele2} eq '-' ) {
        return './.';
    }
    elsif ( $allele1_matches_ref and $allele2_matches_ref ) { # ref and alleles match
        return '0/0'; # homozygous ref
    }
    elsif ( $allele1_matches_ref or $allele2_matches_ref ) { # ref matches one allele, but alleles are not the same
        return '0/1'; # heterozygous ref
    }
    elsif ( $genotype->{allele1} eq $genotype->{allele2} ) { # alleles match, but not the ref
        return '1/1'; # homozygous alt
    }
    else { # not $allele1_matches_ref and not $allele2_matches_ref
        return '1/2';
    }

}

1;

