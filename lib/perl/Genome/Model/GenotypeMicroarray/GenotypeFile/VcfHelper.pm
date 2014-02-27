package Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper { 
    is => 'UR::Singleton',
};

our $supported_info_fields = [
    # SORTED BY ID!
    {
        name => 'allele_frequency',
        id => 'AF',
        header => ',Number=1,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed">',
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
];
sub supported_info_fields {
    return $supported_info_fields;
}

our $info_order;
sub info_order { 
    return $info_order if $info_order;
    $info_order = [ map { $_->{id} } sort { $a->{id} cmp $b->{id} } @$supported_info_fields ];
    return $info_order;
}

our $info_names_and_ids;
sub info_names_and_ids { 
    return $info_names_and_ids if $info_names_and_ids;
    $info_names_and_ids = { map { $_->{name} => $_->{id} } @$supported_info_fields };
    return $info_names_and_ids;
}

our $formats = [qw/ GT OG /];
sub formats {
    return $formats
}

our $format_key_idx = { GT => 0, OG => 1, };
sub format_key_idx {
    return $format_key_idx;
}

sub header {
    my ($self, @sample_names) = @_;

    my @header_lines = (
        '##fileformat=VCFv4.1',
    );
    for my $field ( @$supported_info_fields ) {
        push @header_lines, '##INFO=<ID='.$field->{id}.$field->{header};
    }
    push @header_lines, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
    push @header_lines, '##FORMAT=<ID=AL,Number=1,Type=String,Description="Alleles">';
    push @header_lines, '#'.join("\t", Genome::File::Vcf::Header->column_headers, @sample_names);

    my $header = Genome::File::Vcf::Header->create(lines => \@header_lines);
    return $header;
}

sub convert_genotype_to_vcf_entry {
    my ($self, $genotype) = @_;

    $genotype->{identifiers} = [ $genotype->{id} ];
    $genotype->{chrom} = $genotype->{chromosome};
    $genotype->{reference_allele} = $genotype->{reference};
    $genotype->{alternate_alleles} = $self->alts_for_genotype($genotype);
    $genotype->{quality} = '.';
    $genotype->{_filter} = [];
    $genotype->{info_fields} = $self->info_hash_for_genotype($genotype),
    $genotype->{_format} = $self->formats;
    $genotype->{_format_key_to_idx} = $self->format_key_idx;
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
    for my $info_field ( @{$self->supported_info_fields} ) { 
        my $value = $genotype->{ $info_field->{name} };
        next if not defined $value or $value eq 'NA'; # FIXME convert undef/NA to . ??
        $info_hash{ $info_field->{id} } = $value;
        push @order, $info_field->{id};
    }

    return { hash => \%info_hash, order => \@order, };
}

1;

