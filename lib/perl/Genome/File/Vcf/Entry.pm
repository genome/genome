package Genome::File::Vcf::Entry;

use Carp qw/confess/;
use Genome;

use strict;
use warnings;

# column offsets in a vcf file
use constant {
    CHROM => 0,
    POS => 1,
    ID => 2,
    REF => 3,
    ALT => 4,
    QUAL => 5,
    FILTER => 6,
    INFO => 7,
    FORMAT => 8,
    FIRST_SAMPLE => 9,
};

class Genome::File::Vcf::Entry {
    is => ['UR::Object'],
    has_transient_optional => [
        chrom => {
            is => 'Text',
            doc => 'The sequence name for this site',
        },
        position => {
            is => 'Integer',
            doc => 'The position on the sequence for this site',
        },
        identifiers => {
            is => 'Text',
            doc => 'Optional list Identifiers for this site (e.g., rsid)',
            is_many => 1,
        },
        reference_allele => {
            is => 'Text',
            doc => 'Reference allele at this site',
        },
        alternate_alleles => {
            is => 'Text',
            doc => 'List of alternate alleles at this site',
            is_many => 1,
        },
        quality => {
            is => 'Number',
            doc => 'Quality score for this site',
        },
        filter => {
            is => 'Text',
            doc => 'List of failed filters at this site',
            is_many => 1,
        },
        info_fields => {
            is => 'HASH',
            doc => 'Hash of info fields at this site',
        },
        format => {
            is => 'Text',
            doc => 'Per-sample data format for this site',
            is_many => 1,
        },
        sample_data => {
            is => 'ARRAY',
            doc => 'List of per sample data, with each entry following format',
        },
        # this is built on demand
        _format_key_to_idx => {
            is => 'HASH',
            doc => 'Mapping of format keys to indices, e.g., GT => 0'
        },
        _line => {
            is => 'Text',
        }
    ],
};

sub create {
    my $class = shift;
    return $class->SUPER::create(@_);
}

sub parse {
    my ($self, $line) = @_;
    #my $line = $self->id;
    confess "Attempted to parse null VCF entry" unless $line;

    $self->_line($line);
    my @fields = split("\t", $line);

    # set mandatory fields
    $self->chrom($fields[CHROM]);
    $self->position($fields[POS]);
    $self->identifiers([_parse_list($fields[ID], ',')]);
    $self->reference_allele($fields[REF]);
    $self->alternate_alleles([split(',', $fields[ALT])]);
    $self->quality($fields[QUAL] eq '.' ? undef : $fields[QUAL]);
    $self->filter([_parse_list($fields[FILTER], ',')]);
    $self->info_fields(_parse_info($fields[INFO]));
    $self->format([_parse_list($fields[FORMAT], ':')]);
    $self->sample_data($self->_parse_samples(\@fields)); 
}

# This is to avoid warnings about splitting undef values and to translate
# '.' back to undef.
sub _parse_list {
    my ($identifiers, $delim) = @_;
    return if !$identifiers || $identifiers eq '.';
    return split($delim, $identifiers);
}

# transform a string A=B;C=D;... into a hashref { A=>B, C=>D, ...}
# take care about the string being . (_parse_list does this for us)
sub _parse_info {
    my $info_str = shift;
    my @info_list = _parse_list($info_str, ';');
    return unless @info_list;
    my %info_hash = map {
            my ($k, $v) = split('=', $_, 2);
            $k => $v
        } @info_list;
    return \%info_hash;
}

# this assumes $self->format is set
sub _parse_samples {
    my ($self, $fields) = @_; # arrayref of all vcf fields

    # it is an error to have sample data with no format specification
    confess "VCF entry has sample data but no format specification: " . $self->_line
        if $#$fields >= FIRST_SAMPLE && !defined $self->format;

    my $n_fields = $self->format;

    my @samples;

    for my $i (FIRST_SAMPLE..$#$fields) {
        my $sample_idx = $i - FIRST_SAMPLE + 1;

        my @data = _parse_list($fields->[$i], ':');
        confess "Too many entries in sample $sample_idx:\n"
            ."sample string: $fields->[$i]\n"
            ."format string: " . join(":", $self->format)
            if (@data > $n_fields);

        push(@samples, \@data);
    }
    return \@samples;
}

# Returns a hash of format field names to their indices in per-sample
# lists.
sub _format_field_hash {
    my ($self, $key) = @_;
    my @format = $self->format;
    return unless @format;

    if (!defined $self->_format_key_to_idx) {
        my %h;
        @h{@format} = 0..$#format;
        $self->_format_key_to_idx(\%h);
    }

    return $self->_format_key_to_idx;
}

sub info {
    my ($self, $key) = @_;
    return $self->info_fields unless $key;

    return unless $self->info_fields && exists $self->info_fields->{$key};
    return $self->info_fields->{$key};
}

sub sample_field {
    my ($self, $sample_idx, $field_name) = @_;
    my $sample_data = $self->sample_data;
    confess "Invalid sample index $sample_idx (have $#$sample_data): " . $self->_line
        unless ($sample_idx >= 0 && $sample_idx <= $#$sample_data);

    my $cache = $self->_format_field_hash;
    return unless exists $cache->{$field_name};
    my $field_idx = $cache->{$field_name};
    use Data::Dumper;
    return $sample_data->[$sample_idx]->[$field_idx];
}

sub sample_genotype {
    my ($self, $sample_idx) = @_;
}

sub alleles {
    my $self = shift;
    return ($self->reference_allele, $self->alternate_alleles);
}

sub allelic_distribution {
    my ($self, @sample_indices) = @_;

    # If @sample_indices was not passed, default to all samples
    @sample_indices = 0..$#{$self->sample_data} unless @sample_indices;

    # Find all samples that passed filters
    my @passed_filters = grep {
            my $ft = $self->sample_field($_, "FT");
            !defined $ft || $ft eq "PASS" || $ft eq "."
        } @sample_indices;

    # Get all defined genotypes
    my @gts = grep {defined $_} map {$self->sample_field($_, "GT")} @passed_filters;

    # Split gts on |/ and flatten into one list
    my @allele_indices = grep {defined $_ && $_ ne '.'} map {split("[/|]", $_)} @gts;

    # Get list of all alleles (ref, alt1, ..., altn)
    my @alleles = $self->alleles;
    my %counts;
    my $total = 0;
    for my $a (@alleles[@allele_indices]) {
        ++$counts{$a};
        ++$total;
    }
    return ($total, %counts);
}

1;
