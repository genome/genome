package Genome::File::Vcf::VepConsequenceParser;

# This module parses the CSQ info field added by the VEP annotator using the
# format specification in the vcf header.

use Data::Dumper;
use Carp qw/confess/;

use strict;
use warnings;

# Constructor
#  Args:
#   vcf_header      - a Genome::File::Vcf::Header
#   options         - a hash with optional arguments:
#       tag_name    - the name of the vep consequence tag (DEFAULT: CSQ)
sub new {
    my ($class, $vcf_header, %options) = @_;

    my $filters = delete $options{filters} || [];
    my $tag_name = delete $options{tag_name} || "CSQ";
    confess "Unknown options: " . Dumper(\%options) . " to " . __PACKAGE__ . " constructor"
        if (%options);

    my $csq_info = $vcf_header->info_types->{$tag_name};
    confess "INFO tag $tag_name not found in VCF header" unless $csq_info;

    my $format_str = $csq_info->{description};
    $format_str =~ s/^[^|]* Format: //g;

    confess "Failed to extract format string from info description for $tag_name: " . Dumper($csq_info)
        unless $format_str;

    my @format = map {lc} split("\\|", $format_str);

    my $self = {
        filters => $filters,
        tag_name => $tag_name,
        fields => \@format,
        field_index => [map {$format[$_] => $_} 0..$#format],
    };

    bless $self, $class;
    return $self;
}

# In the presence of indels, VEP just chops off the first base of all alts
# and refs, which isn't precisely the right thing to do. We can't just prepend
# the reference base to all alts, we need to construct a map from what they
# will call the allele to what it is in the original input file.
sub resolve_alleles {
    my ($entry) = @_;
    if ($entry->has_indel) {
        return map {substr($_, 1) || '-' => $_} @{$entry->{alternate_alleles}};
    }
    return;
}

sub _passes_filters {
    my ($self, $annotation) = @_;
    for my $f (@{$self->{filters}}) {
        if (!$f->($annotation)) {
            return;
        }
    }
    return 1;
}

sub process_entry {
    my ($self, $entry) = @_;

    my $value = $entry->info($self->{tag_name});
    return unless $value;

    my @annotations = split(",", $value);

    my %allele_map = resolve_alleles($entry);

    my $rv = {};
    for my $ann (@annotations) {
        my @values = split("\\|", $ann, -1);
        my %h = map {$self->{fields}[$_] => $values[$_]} 0..$#values;
        my $allele;
        if (%allele_map) {
            if (!exists $allele_map{$h{allele}}) {
                confess "Unknown allele from vep in vcf entry: $h{allele} not found in entry "
                    . Dumper($entry);
            }

            $allele = $allele_map{$h{allele}};
        }
        else {
            $allele = $h{allele};
        }
        if ($self->_passes_filters(\%h)) {
            push @{$rv->{$allele}}, \%h;
        }
    }

    return $rv;
}

sub _format_transcript {
    my ($self, $transcript) = @_;
    return join("|", map {$transcript->{$_} || ''} @{$self->{fields}});
}

sub format_transcripts {
    my ($self, $transcripts) = @_;
    return '.' unless %$transcripts;
    my @flat_transcripts = map { @$_ } values %$transcripts;
    return join(",", map {$self->_format_transcript($_)} @flat_transcripts);
}

sub transcripts {
    my ($self, $entry, $allele) = @_;

    my $processed_entry = $self->process_entry($entry);
    my $transcripts = $processed_entry->{$allele};
    return @{$transcripts};
}

sub canonical_transcripts {
    my ($self, $entry, $allele) = @_;

    my @transcripts = @{$self->process_entry($entry)->{$allele}};
    return grep {$_->{'canonical'} eq 'YES'} @transcripts;
}

1;
