package Genome::File::Vcf::DbsnpAFParser;

use strict;
use warnings;
use Genome;

use constant TAG => "CAF";

sub new {
    my ($class, $vcf_header) = @_;
    my $tag = TAG;
    my $caf_info = $vcf_header->info_types->{$tag};
    die sprintf("INFO tag %s not found in VCF header", TAG) unless $caf_info;
    my $self = {};

    bless $self, $class;
    return $self;
}

sub process_entry {
    my ($self, $entry) = @_;

    my $caf = $entry->info(TAG);
    unless (defined $caf) {
        return;
    }
    my ($caf_string) = $caf =~ /\[([0-9,\.]*)\]/;
    unless ($caf_string) {
        die sprintf("Invalid CAF entry '%s'\n", $caf);
    }

    my @allele_frequencies = split ",", $caf_string;
    unless (scalar @allele_frequencies == scalar ($entry->alleles)) {
        die sprintf("Frequency list and allele list differ in length. %s:(%s) %s:(%s)",
            scalar @allele_frequencies, join(",", @allele_frequencies),
            scalar ($entry->alleles), join(",", $entry->alleles));
    }

    my %af;
    @af{$entry->alleles} = @allele_frequencies;

    return \%af;
}

1;

