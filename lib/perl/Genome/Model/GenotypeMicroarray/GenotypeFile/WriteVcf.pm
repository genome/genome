package Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Header;
use Genome::File::Vcf::Entry;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf { 
    is => 'Genome::Utility::IO::Writer',
};

my @genotype_info_fields = (
    # SORTED BY ID!
    {
        name => 'allele_count',
        id => 'AC', 
        header => ',Number=1,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">', 
    },
    {
        name => 'allele_frequency',
        id => 'AF',
        header => ',Number=1,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed',
    },
    {
        name => 'allele_number',
        id => 'AN',
        header => ',Number=1,Type=Integer,Description="Total number of alleles in called genotypes',
    },
    {
        name => 'cnv_confidence',
        id => 'CC',
        header => ',Number=1,Type=Float,Description="CNV Confidence">',
    },
    {
        name => 'cnv_value',
        id => 'CV',
        header => ',Number=1,Type=Float,Description="CNV Value">',
    },
    {
        name => 'log_r_ratio',
        id => 'LR',
        header => ',Number=1,Type=Float,Description="Log R Ratio">',
    },
    {
        name => 'gc_score',
        id => 'GC',
        header => ',Number=1,Type=Float,Description="GC Score">',
    },
    {
        name => 'alleles', # origninal_genotype
        id => 'OG',
        header => ',Number=1,Type=Integer,Description="Original genotype calls">', 
    },
);


sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $header = $self->_construct_header_string;
    return if not $header;

    $self->output->print($header);

    return $self;
}

sub write_one {
    my ($self, $genotype) = @_;

    Carp::confess('No genotype to write!') if not $genotype;

    $self->output->print( 
        join(
            "\t", 
            $genotype->{chromosome},                        # CHROM
            $genotype->{position},                          # POS
            $genotype->{id},                                # ID
            $genotype->{reference},                         # REF
            $self->alt_string_for_genotype($genotype),      # ALT
            '.',                                            # QUAL
            '.',                                            # FILTER
            $self->info_string_for_genotype($genotype),     # INFO
            $self->format_string_for_genotype($genotype),   # FORMAT
        )."\n"
    );

    return 1;
}

sub _construct_header_string {
    my $self = shift;

    # META
    my $header = "##fileformat=VCFv4.1\n";

    # INFO
    for my $field ( sort { $a->{id} cmp $b->{id} } @genotype_info_fields ) {
        $header .= '##INFO=<ID='.$field->{id}.$field->{header}."\n";
    }

    # FILTER...none
    # FORMAT...none

    # HEADERS
    $header .= '#'.join("\t", (qw/ CHROM POS ID REF ALT QUAL FILTER INFO FORMAT /))."\n"; # FIXME SAMPLENAME

    return $header;
}

sub alt_string_for_genotype {
    my ($self, $genotype) = @_;

    Carp::confess('No genotype to get alt string!') if not $genotype;
    # check allele1 allele2 ref?

    # alleles are dashes
    return '.' if $genotype->{allele1} eq '-' and $genotype->{allele2} eq '-';

    # see which match the ref
    my @alts;
    for my $allele_key ( map{ 'allele'.$_ } (1..2) ) {
        next if $genotype->{$allele_key} eq '-';
        push @alts, $genotype->{$allele_key} if $genotype->{$allele_key} ne $genotype->{reference};
    }

    # both match!
    return '.' if not @alts;
    # retrun alts
    return join(',', sort { $a cmp $b } @alts);
}

sub info_string_for_genotype {
    my ($self, $genotype) = @_;

    Carp::confess('No genotype to get info string!') if not $genotype;

    # key=value[;key=value]
    my @infos;
    for my $field ( sort { $a->{id} cmp $b->{id} } @genotype_info_fields ) {
        next if not defined $genotype->{ $field->{name} } or $genotype->{ $field->{name} } eq 'NA';
        push @infos, $field->{id}.'='.$genotype->{ $field->{name} };
    }

    return '.' if not @infos;
    return join(';', @infos);
}

sub format_string_for_genotype {
    my ($self, $genotype) = @_;

    Carp::confess('No genotype to get info string!') if not $genotype;

    my $allele1_matches_ref = $genotype->{allele1} eq $genotype->{reference};
    my $allele2_matches_ref = $genotype->{allele2} eq $genotype->{reference};

    my @formats = 'GT';
    if ( $genotype->{allele1} eq '-' and $genotype->{allele2} eq '-' ) {
        push @formats, './.';
    }
    elsif ( $allele1_matches_ref and $allele2_matches_ref ) { # ref and alleles match
        push @formats, '0/0'; # homozygous ref
    }
    elsif ( $allele1_matches_ref or $allele2_matches_ref ) { # ref matches one allele, but alleles are not the same
        push @formats, '0/1'; # heterozygous ref
    }
    elsif ( $genotype->{allele1} eq $genotype->{allele2} ) { # alleles match, but not the ref
        push @formats, '1/1'; # homozygous alt
    }
    else { # not $allele1_matches_ref and not $allele2_matches_ref
        push @formats, '1/2';
    }

    return join("\t", @formats);
}

1;

