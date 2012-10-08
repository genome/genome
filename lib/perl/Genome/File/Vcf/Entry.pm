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

sub new {
    my ($class, $header, $line) = @_;
    my $self = {
        header => $header,
        _line => $line,
    };
    bless $self, $class;
    $self->_parse;
    return $self;
}

sub _parse {
    my ($self) = @_;
    confess "Attempted to parse null VCF entry" unless $self->{_line};

    my @fields = split("\t", $self->{_line});
    $self->{_fields} = \@fields;

    # set mandatory fields
    $self->{chrom} = $fields[CHROM];
    $self->{position} = $fields[POS];

    $self->{identifiers} = undef;
    my @identifiers = _parse_list($fields[ID], ',');
    $self->{identifiers} = \@identifiers if @identifiers;

    $self->{reference_allele} = $fields[REF];
    $self->{alternate_alleles} = [split(',', $fields[ALT])];
    $self->{quality} = $fields[QUAL] eq '.' ? undef : $fields[QUAL];
    $self->{filter} = [_parse_list($fields[FILTER], ',')];
    $self->{format} = [_parse_list($fields[FORMAT], ':')];
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
        # we do this rather than just split to handle Flag types (they have undef values)
        $k => $v
        } @info_list;
    return \%info_hash;
}

# this assumes $self->{format} is set
sub _parse_samples {
    my ($self) = @_; # arrayref of all vcf fields

    my $fields = $self->{_fields};

    # it is an error to have sample data with no format specification
    confess "VCF entry has sample data but no format specification: " . $self->{_line}
        if $#$fields >= FIRST_SAMPLE && !defined $self->{format};

    my $n_fields = $self->{format};

    my @samples;

    for my $i (FIRST_SAMPLE..$#$fields) {
        my $sample_idx = $i - FIRST_SAMPLE + 1;

        my @data = _parse_list($fields->[$i], ':');
        confess "Too many entries in sample $sample_idx:\n"
            ."sample string: $fields->[$i]\n"
            ."format string: " . join(":", $self->{format})
            if (@data > $n_fields);

        push(@samples, \@data);
    }
    return \@samples;
}

# Returns a hash of format field names to their indices in per-sample
# lists.
sub format_field_index {
    my ($self, $key) = @_;
    if (!exists $self->{_format_key_to_idx}) {
        my @format = @{$self->{format}};
        return unless @format;

        my %h;
        @h{@format} = 0..$#format;
        $self->{_format_key_to_idx} = \%h;
    }

    return $self->{_format_key_to_idx};
}

sub info {
    my ($self, $key) = @_;
    $self->{info_fields} = _parse_info($fields[INFO]) unless exists $self->{info_fields};

    return $self->{info_fields} unless $key;

    return unless $self->{info_fields} && exists $self->{info_fields}->{$key};

    # The 2nd condition is to deal with flags, they may exist but not have a value
    return $self->{info_fields}->{$key} || exists $self->{info_fields}->{$key};
}

sub info_for_allele {
    my ($self, $allele, $key) = @_;

    # nothing to return
    return unless defined $self->{info_fields};

    # we don't have that allele, or it is the reference (idx 0)
    my $idx = $self->allele_index($allele);
    return unless defined $idx && $idx > 0;
    --$idx; # we don't care about the reference allele

    # no header! what are you doing?
    confess "info_for_allele called on entry with no vcf header!" unless $self->{header};
    
    my @keys = defined $key ? $key : keys %{$self->{info_fields}};

    my %result;
    for my $k (@keys) {
        my $type = $self->{header}->info_types->{$k};
        warn "Unknown info type $k encountered!" if !defined $type;
        if (defined $type && $type->{number} eq 'A') { # per alt field

            my @values = split(',', $self->{info_fields}->{$k});
            next if $idx > $#values;
            $result{$k} = $values[$idx];
        } else {
            $result{$k} = $self->{info_fields}->{$k}
        }
    }
    my $rv = $key ? $result{$key} : \%result;
    return $rv;
}

sub sample_data {
    my $self = shift;
    if (!exists $self->{_sample_data}) {
        $self->{_sample_data} = $self->_parse_samples;
    }
    return $self->{_sample_data};
}

sub sample_field {
    my ($self, $sample_idx, $field_name) = @_;
    my $sample_data = $self->sample_data;
    confess "Invalid sample index $sample_idx (have $#$sample_data): " . $self->{_line}
        unless ($sample_idx >= 0 && $sample_idx <= $#$sample_data);

    my $cache = $self->format_field_index;
    return unless exists $cache->{$field_name};
    my $field_idx = $cache->{$field_name};
    return $sample_data->[$sample_idx]->[$field_idx];
}

sub alleles {
    my $self = shift;
    return ($self->{reference_allele}, @{$self->{alternate_alleles}});
}

# return the index of the given allele, or undef if not found
# note that 0 => reference, 1 => first alt, ...
sub allele_index {
    my ($self, $allele) = @_;
    my @a = $self->alleles;
    my @idx = grep {$a[$_] eq $allele} 0..$#a;
    return unless @idx;
    return $idx[0]
}

sub allelic_distribution {
    my ($self, @sample_indices) = @_;

    my $format = $self->format_field_index;
    my $gtidx = $format->{GT};
    my $sample_data = $self->sample_data;
    return unless defined $gtidx && defined $sample_data;

    # If @sample_indices was not passed, default to all samples
    @sample_indices = 0..$#{$sample_data} unless @sample_indices;

    # Find all samples that passed filters
    if (exists $format->{FT}) {
        my $ftidx = $format->{FT};
        @sample_indices = grep {
            my $ft = $sample_data->[$_]->[$ftidx];
            !defined $ft || $ft eq "PASS" || $ft eq "."
            } @sample_indices;
    }

    # Get all defined genotypes
    my @gts = 
        grep {defined $_} 
        map {$sample_data->[$_]->[$gtidx]}
        @sample_indices;

    # Split gts on |/ and flatten into one list
    my @allele_indices = 
        grep {defined $_ && $_ ne '.'}
        map {split("[/|]", $_)}
        @gts;

    # Get list of all alleles (ref, alt1, ..., altn)
    my %counts;
    my $total = 0;
    for my $a (@allele_indices) {
        ++$counts{$a};
        ++$total;
    }
    return ($total, %counts);
}

sub is_filtered {
    my $self = shift;
    return $self->{filter} && grep { $_ && $_ ne "PASS" && $_ ne "."} @{$self->{filter}};
}

1;
