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


    my %info_type_numbers_and_subclasses = (
        '.' => 'Dot',
        'R' => 'R',
    );
    my $info_type_number = $caf_info->{number};
    die 'Invalid CAF info type number! '.$info_type_number if not exists $info_type_numbers_and_subclasses{$info_type_number};
    
    my $sub_class = $class.'ForInfoType'.$info_type_numbers_and_subclasses{$info_type_number};

    return bless {}, $sub_class;
}

package Genome::File::Vcf::DbsnpAFParserForInfoTypeDot;

use parent 'Genome::File::Vcf::DbsnpAFParser';

use constant TAG => "CAF";

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

package Genome::File::Vcf::DbsnpAFParserForInfoTypeR;

use parent 'Genome::File::Vcf::DbsnpAFParser';

sub process_entry {
    my ($self, $entry) = @_;
    my %af;
    for my $allele ( $entry->alleles ) {
        my $af = $entry->info_for_allele($allele, "CAF");
        $af{$allele} = $af;
    }

    return \%af;
}

1;

