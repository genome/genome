package Genome::File::Vcf::Entry;

use Data::Dumper;
use Carp qw/confess/;
use Genome;
use List::Util qw/first/;

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

=head1 NAME

Genome::File::Vcf::Entry - Representation of a single (non-header) line in a
VCF file.

=head1 SYNOPSIS

    use Genome::File::Vcf::Entry;
    use Genome::File::Vcf::Reader;

    my $reader = new Genome::File::Vcf::Reader("input.vcf");
    my $entry = $reader->next;

    printf "entry: %s\n", $entry->to_string;

    printf "CHROM:   %s\n", $entry->{chrom};
    printf "POS:     %s\n", $entry->{position};
    printf "IDENT:   %s\n", join(",", $entry->{identifiers});
    printf "REF:     %s\n", $entry->{reference_allele};
    printf "ALT:     %s\n", join(",", $entry->{alternate_alleles});
    printf "QUAL:    %s\n", $entry->{quality};
    # note that the following involve method calls rather than properties
    printf "FILTERS: %s\n" . join(";", $entry->filters);
    ...

    my $gmaf = $entry->info("GMAF");
    if (defined $gmaf) {
        ...
    }

    # REF: A, ALT: AC, ACT
    # Get all info fields for allele AC, paying attention to per-alt fields.
    my $ac_info_all = $entry->info_for_allele("AC");

    # or just get one field
    my $ac_info_gmaf = $entry->info_for_allele("AC", "GMAF");


    # Assume the format fields for the entry ar GT:DP:FT
    my @format_names = $entry->format;
    # @format_names = ("GT", "DP", "FT")

    my $idx = $entry->format_field_index("DP"); # $idx = 1
    my $idx = $entry->format_field_index("GL"); # $idx = undef

    my $cache = $entry->format_field_index;
    # $cache = {"GT" => 0, "DP" => 1, "FT" => 2}

=head1 DESCRIPTION

=head2 Public data members

=over 12

=item C<< $entry->{chrom} >>

    chromosome or sequence name

=item C<< $entry->{position} >>

    position on the chromosome

=item C<< $entry->{identifiers} >>

    arrayref of identifiers (e.g., rsids) for this entry

=item C<< $entry->{reference_allele} >>

    the reference allele for this entry

=item C<< $entry->{alternate_alleles} >>

    arrayref of alternate alleles for this entry

=item C<< $entry->{quality} >>

    phred scaled quality score

=head2 Methods

=back

=over 12

=item C<new>

Constructs a new Genome::File::Vcf::Entry from the given header and raw text
line.

=cut

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

    my @identifiers = _parse_list($fields[ID], ',');
    $self->{identifiers} = \@identifiers;

    $self->{info_fields} = undef; # we parse them lazily later

    $self->{reference_allele} = $fields[REF];
    my @alts = _parse_list($fields[ALT], ',');
    $self->{alternate_alleles} = \@alts;
    $self->{quality} = $fields[QUAL] eq '.' ? undef : $fields[QUAL];
    $self->{_filter} = [_parse_list($fields[FILTER], ',')];
    $self->{_format} = [_parse_list($fields[FORMAT], ':')];
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
    my %info_hash;
    my @info_order;
    for my $i (@info_list) {
        my ($k, $v) = split('=', $i, 2);
        # we do this rather than just split to handle Flag types (they have undef values)
        $info_hash{$k} = $v;
        push(@info_order, $k);
    }
    return { hash => \%info_hash, order => \@info_order };
}

# this assumes $self->{_format} is set
sub _parse_samples {
    my ($self) = @_; # arrayref of all vcf fields

    my $fields = $self->{_fields};

    # it is an error to have sample data with no format specification
    confess "VCF entry has sample data but no format specification: " . $self->{_line}
        if $#$fields >= FIRST_SAMPLE && !defined $self->{_format};

    my $n_fields = $#{$self->{_format}} + 1;

    my @samples;

    for my $i (FIRST_SAMPLE..$#$fields) {
        my $sample_idx = $i - FIRST_SAMPLE + 1;

        my @data = _parse_list($fields->[$i], ':');
        confess "Too many entries in sample $sample_idx:\n"
            ."sample string: $fields->[$i]\n"
            ."format string: " . join(":", $self->{_format})
            if (@data > $n_fields);

        push(@samples, \@data);
    }
    return \@samples;
}

=item C<alleles>

Returns a list of all alleles (reference + alternates) for this entry. This
list can be used to convert GT sample indexes to actual allele values.

=cut

sub alleles {
    my $self = shift;
    return ($self->{reference_allele}, @{$self->{alternate_alleles}});
}

=item C<has_indel>

Returns true if the given entry has an insertion or deletion, false otherwise.

=cut

sub has_indel {
    my $self = shift;
    for my $alt (@{$self->{alternate_alleles}}) {
        if (length($alt) != length($self->{reference_allele})) {
            return 1;
        }
    }
    return 0;
}

=item C<allele_index>

Returns the index of the given allele, or undef if not found.
Note that 0 => reference, 1 => first alt, ...

params:
    $allele - the allele to find

=cut

sub allele_index {
    my ($self, $allele) = @_;
    my @a = $self->alleles;
    my @idx = grep {$a[$_] eq $allele} 0..$#a;
    return unless @idx;
    return $idx[0]
}

=item C<allelic_distribution>

Examines GT fields in sample data and returns (#total alleles, counts)
where counts is a hash mapping allele index -> frequency.

=cut

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

=item C<info>

With no arguments, returns a hashref of all info fields for this entry.
With a single argument, returns just the given info field.

params:
    $key - optional info field name (all are returned if not given)

See info_for_allele for an example.

=cut

sub info {
    my ($self, $key) = @_;

    if (!$self->{info_fields}) {
        $self->{info_fields} = _parse_info($self->{_fields}->[INFO]);
    }

    my $hash = $self->{info_fields}{hash};
    return $hash unless $key;

    return unless $hash && exists $hash->{$key};

    # The 2nd condition is to deal with flags, they may exist but not have a value
    return $hash->{$key} || exists $hash->{$key};
}

=item C<info_for_allele>

Similar to the info method, but if there are any I<per-allele> info fields,
only the values that apply to the allele in question are included.

params:
    $allele - the alternate allele of interest
    $key - optional info field name (all are returned if not given)

=head3 Example:

    # X is per-allele, Y is not.
    ... REF ALT  ... INFO
    ... AC  A,AT ... X=2,3;Y=4

    $entry->info == (X => "2,3", Y => 4)
    $entry->info("X") == "2,3"

    $entry->info_for_allele("A") == (X => 2, Y => 4)
    $entry->info_for_allele("A", "X") == 2

    $entry->info_for_allele("AT") == (X => 3, Y => 4)
    $entry->info_for_allele("AT", "X") == 3

=cut

sub info_for_allele {
    my ($self, $allele, $key) = @_;

    # nothing to return
    my $hash = $self->info;
    return unless defined $hash;

    # we don't have that allele, or it is the reference (idx 0)
    my $idx = $self->allele_index($allele);
    return unless defined $idx && $idx > 0;
    --$idx; # we don't care about the reference allele

    # no header! what are you doing?
    confess "info_for_allele called on entry with no vcf header!" unless $self->{header};

    my @keys = defined $key ? $key : keys %$hash;

    my %result;
    for my $k (@keys) {
        my $type = $self->{header}->info_types->{$k};
        warn "Unknown info type $k encountered!" if !defined $type;
        if (defined $type && $type->{number} eq 'A') { # per alt field

            my @values = split(',', $hash->{$k});
            next if $idx > $#values;
            $result{$k} = $values[$idx];
        } else {
            $result{$k} = $hash->{$k}
        }
    }
    my $rv = $key ? $result{$key} : \%result;
    return $rv;
}

# remove PASS if the entry is filtered
sub _check_filters {
    my $self = shift;
    if ($self->is_filtered) {
        @{$self->{_filter}} = grep {!/^PASS$/} @{$self->{_filter}};
    }
}

=item C<filters>

Return a list of filters that this site has failed.

=cut

sub filters {
    my ($self) = @_;
    return @{$self->{_filter}};
}

=item C<clear_filters>

Remove any site filters from the entry.

=cut

sub clear_filters {
    my $self = shift;
    $self->{_filter} = [];
}

=item C<add_filter>

Add the given filter to the entry (indicates that the entire site is filtered).

params:
    $filter - the name of the filter to apply

=cut

sub add_filter {
    my ($self, $filter) = @_;
    push @{$self->{_filter}}, $filter;
    $self->_check_filters;
}

=item C<is_filtered>

Returns true if the entire site represented by this entry is filtered, indicated
by anything other than PASS or . in the FILTER column.

=cut

sub is_filtered {
    my $self = shift;
    return grep { $_ && $_ ne "PASS" && $_ ne "."} @{$self->{_filter}};
}

sub _prepend_format_field {
    my ($self, $field) = @_;

    return 0;
}

=item C<add_format_field>

Add the given field to the FORMAT specification for this entry and return its
index. If the field previously existed, its index is returned. Typically, fields
are appended to the FORMAT list but "GT" is a special case: it is prepended as
per the vcf spec. This method shifts all sample data appropriately when this
happens.

params:
    $field - the field to add. This B<should> exist in the VCF header.

=head3 Example:

    # FORMAT = GL
    $entry->add_format_field("DP") # returns 1
    # FORMAT = GL:DP
    $entry->add_format_field("GT") # returns 0
    # FORMAT = GT:GL:DP

=cut

sub add_format_field {
    my ($self, $field) = @_;
    if (!exists $self->{header}->format_types->{$field}) {
        confess "Format field '$field' does not exist in header";
    }

    if (exists $self->{_format_key_to_idx}{$field}) {
        return $self->{_format_key_to_idx}{$field};
    }
    elsif ($field eq "GT") {
        # GT is a special case in the Vcf spec. If present, it must be the first field

        # Increment all indexes by one since we are prepending something
        my $fmtidx = $self->{_format_key_to_idx};
        %$fmtidx = map {$_ => $fmtidx->{$_} + 1} keys %$fmtidx;

        # Prepend undef to all the existing sample entries
        unshift @{$self->{_format}}, $field;
        my $sample_data = $self->sample_data;
        for my $i (0..$#$sample_data) {
            unshift @{$sample_data->[$i]}, undef;
        }

        $self->{_format_key_to_idx}{$field} = 0;
        return 0;
    }
    else {
        my $idx = scalar @{$self->{_format}};
        push @{$self->{_format}}, $field;
        $self->{_format_key_to_idx}{$field} = $idx;
        return $idx;
    }
}

=item C<format>

Returns the list of FORMAT fields for this entry.

=cut

sub format {
    my $self = shift;
    return @{$self->{_format}};
}

=item C<format_field_index>

If called with no parameters, returns a hash of format field names to their
indexes in the FORMAT specification for this entry.

If called with a single parameter, it is interpreted as a format field name
and its index is returned (or undef if it does not exist). See L<SYNOPSIS> for
an example.

=cut

sub format_field_index {
    my ($self, $key) = @_;
    if (!exists $self->{_format_key_to_idx}) {
        my @format = @{$self->{_format}};
        return unless @format;

        my %h;
        @h{@format} = 0..$#format;
        $self->{_format_key_to_idx} = \%h;
    }

    return $self->{_format_key_to_idx}{$key} if $key;
    return $self->{_format_key_to_idx};
}

=item C<sample_data>

Get a raw reference to all sample data for this entry. Sample data is lazily
parsed, so if this method is never called, you never pay for the parsing. Making
modifications to the return value can corrupt the internal state of the Entry
object. Be careful.

returns: $x = [ [..sample0 values..], ..., [..sampleN values..] ]

=cut

sub sample_data {
    my $self = shift;
    if (!exists $self->{_sample_data}) {
        $self->{_sample_data} = $self->_parse_samples;
    }
    return $self->{_sample_data};
}

=item C<sample_field>

Get sample data values.

params:
    $sample_idx - sample index
    $field_name - name of the value to fetch (e.g., GT, DP)

=head3 Example:
    FORMAT    SAMPLE1      SAMPLE2
    GT:DP:FT  0/1:23:PASS  1/1:24:FalsePositive

    $entry->sample_field(0, "GT") == "0/1"
    $entry->sample_field(0, "DP") == 23
    $entry->sample_field(1, "FT") == "FalsePositive"

=cut

sub sample_field {
    my ($self, $sample_idx, $field_name) = @_;
    my $n_samples = scalar $self->{header}->sample_names;
    confess "Invalid sample index $sample_idx (have $n_samples samples): "
        unless ($sample_idx >= 0 && $sample_idx <= $n_samples);

    my $sample_data = $self->sample_data;
    return unless $sample_idx <= $#$sample_data;

    my $cache = $self->format_field_index;
    return unless exists $cache->{$field_name};
    my $field_idx = $cache->{$field_name};
    return $sample_data->[$sample_idx]->[$field_idx];
}

sub _extend_sample_data {
    my ($self, $idx) = @_;

    my $sample_data = $self->sample_data;
    return if $idx <= $#$sample_data;

    # We might have some undefined entries between here and the target sample
    # index. We should fill those with . (or ./. for GT).
    my $n_fields = scalar @{$self->{_format}};

    my @empty_data = (undef) x $n_fields;
    for my $i ($#$sample_data+1..$idx) {
        push @$sample_data, [@empty_data];
    }
}

=item C<set_sample_field>

Sets the value of per-sample fields. It is assumed that the $field_name already
exists in the FORMAT specification for this entry. If this is not the case, it
can be added with $entry->add_format_field.

params:
    $sample_idx - sample index
    $field_name - name of the field to set (must exist in FORMAT)
    $value      - the value to set $field_name to

=cut

sub set_sample_field {
    my ($self, $sample_idx, $field_name, $value) = @_;

    my $sample_data = $self->sample_data;
    my $n_samples = scalar $self->{header}->sample_names;
    confess "Invalid sample index $sample_idx (have $n_samples samples): "
        unless ($sample_idx >= 0 && $sample_idx <= $n_samples);

    my $cache = $self->format_field_index;
    if (!exists $cache->{$field_name}) {
        confess "Unknown format field $field_name";
    }

    $self->_extend_sample_data($sample_idx);
    my $field_idx = $cache->{$field_name};
    $sample_data->[$sample_idx]->[$field_idx] = $value;
}

=item C<filter_calls_involving_only>

Applies a per-sample filter (FT tag) to any genotype calls that involve only the
specified alleles.

params (pass by hash):
    filter_name - the name of the filter to apply
    alleles - arrayref of alleles

=head3 Example:

BEFORE:
    FORMAT    SAMPLE1      SAMPLE2
    GT:DP:FT  0/1:23:PASS  0/0:24:.

$entry->filter_calls_involving_only(
    filter_name => "WILDTYPE",
    alleles => [$entry->{reference_allele}]
    );

AFTER:
    BEFORE:
    FORMAT    SAMPLE1      SAMPLE2
    GT:DP:FT  0/1:23:PASS  0/0:24:WILDTYPE

Sample 2 is filtered because its call involves only the reference allele.

=cut

sub filter_calls_involving_only {
    my ($self, %options) = @_;
    my $filter_name = delete $options{filter_name} || confess "Missing argument: filter_name";
    my $alleles = delete $options{alleles} || confess "Missing argument: alleles";

    my $gt_idx = $self->format_field_index("GT");
    # No genotypes here, just bail.
    return unless defined $gt_idx;

    my $sample_data = $self->sample_data;
    # No sample data
    return unless $#$sample_data >= 0;

    # This will return the existing index if the thing already exists
    my $ft_idx = $self->add_format_field("FT");


    my @alleles = $self->alleles;
    my %filter_alleles = map {$_ => undef} @$alleles;

    for my $sample_idx (0..$#$sample_data) {
        my $gt = $self->sample_field($sample_idx, "GT");
        next unless defined $gt;

        my @alt_indices = split("[/|]", $gt);
        my @call_alleles = map {$alleles[$_]} @alt_indices;
        # We filter the call when there is nothing in @call_alleles that is not
        # in %filter_alleles
        if (!defined first { !exists $filter_alleles{$_} } @call_alleles) {
            # This will return the existing index if the thing already exists
            $self->set_sample_field($sample_idx, "FT", $filter_name);
        }
    }
}

sub genotype_for_sample {
    my ($self, $sample_index) = @_;
    my @vars = $self->{alternate_alleles};
    my $gt = $self->sample_field($sample_index, 'GT');
    unless ($gt) {
        confess "No sample for index $sample_index";
    }
    return Genome::File::Vcf::Genotype->new($self->{reference_allele}, \@vars, $gt);
}

=item C<to_string>

Returns a string representation of the entry in VCF format.

=back

=cut


sub to_string {
    my ($self) = @_;

    my %info = %{$self->info};
    my %info_keys = map {$_ => undef} keys %info;

    my @info_order;
    for my $info_key (@{$self->{info_fields}{order}}) {
        push(@info_order, $info_key) if exists $info{$info_key};
        delete $info_keys{$info_key};
    }

    push(@info_order, keys %info_keys);

    # We want to display ./. when the GT format field is undefined,
    # but . otherwise. To this end, we build an array containing the
    # proper string to represent undef for each info field present.
    my @format_undef = ('.') x scalar @{$self->{_format}};
    my $gtidx = $self->format_field_index("GT");
    if (defined $gtidx) {
        $format_undef[$gtidx] = './.';
    }

    return join("\t",
        $self->{chrom},
        $self->{position},
        join(",", @{$self->{identifiers}}) || '.',

        $self->{reference_allele} || '.',
        join(",", @{$self->{alternate_alleles}}) || '.',
        $self->{quality} || '.',
        join(";", @{$self->{_filter}}) || '.',
        join(";",
            map {
                defined $info{$_} ?
                    join("=", $_, $info{$_})
                    : $_
            } @info_order) || '.',
        join(":", @{$self->{_format}}) || '.',
        map {
            # Join values for an individual sample
            my $values = $_;
            join(":", map {
                defined $values->[$_] ? $values->[$_] : $format_undef[$_]
                } 0..$#$values) || '.'
            } @{$self->sample_data}, # for each sample
        );
}



1;
