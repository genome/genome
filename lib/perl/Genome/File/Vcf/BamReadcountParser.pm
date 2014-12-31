package Genome::File::Vcf::BamReadcountParser;

use strict;
use warnings;

use Genome::File::BamReadcount::Entry;
use JSON qw(to_json from_json);
use Scalar::Util qw(looks_like_number);


sub get_bam_readcount_entries {
    my ($vcf_entry, $sample_idx) = @_;

    my $string = $vcf_entry->sample_field($sample_idx, 'BRCT');
    unless ($string) {
        return {};
    }
    my $per_allele_hash = from_json(decode($string));

    my %readcount_entries;
    for my $offset (grep {looks_like_number($_)} keys %$per_allele_hash) {
        my $bam_readcount_string = $per_allele_hash->{$offset};
        die sprintf("No entry for offset $offset in entry %s", $vcf_entry->to_string) unless $bam_readcount_string;
        if ($bam_readcount_string eq '.') {
            $readcount_entries{$offset} = undef;
        }
        else {
            $readcount_entries{$offset} = Genome::File::BamReadcount::Entry->new($bam_readcount_string);
        }
    }

    my %readcount_entries_per_allele;
    for my $alt_allele (@{$vcf_entry->{alternate_alleles}}) {
        my $offset = $per_allele_hash->{$alt_allele};
        die sprintf("No BRCT string for allele %s in entry %s", $alt_allele, $vcf_entry->to_string) unless defined $offset;
        $readcount_entries_per_allele{$alt_allele} = $readcount_entries{$offset};
    }

    return \%readcount_entries_per_allele;
}

sub add_bam_readcount_entries {
    my ($vcf_entry, $reader_hash, $readcount_tag) = @_;

    unless (@{$vcf_entry->{alternate_alleles}}) {
        return;
    }
    $vcf_entry->add_format_field($readcount_tag);
    my $allele_offsets = get_allele_offsets($vcf_entry);
    for my $sample_idx (keys %$reader_hash) {
        my $readcount_string = generate_readcount_string(
            $vcf_entry,
            $allele_offsets,
            $reader_hash->{$sample_idx}
        );
        $vcf_entry->set_sample_field($sample_idx, $readcount_tag, $readcount_string);
    }
}

sub get_allele_offsets {
    my $vcf_entry = shift;
    my $offsets;
    for my $alt_allele (@{$vcf_entry->{alternate_alleles}}) {
        if ($vcf_entry->is_deletion($alt_allele)) {
            $offsets->{$alt_allele} = 1;
        }
        else {
            $offsets->{$alt_allele} = 0;
        }
    }
    return $offsets;
}

sub generate_readcount_string {
    my ($vcf_entry, $allele_offsets, $reader) = @_;
    my %allele_hash = %$allele_offsets;
    my $offsets = Set::Scalar->new(values %allele_hash);
    for my $offset ($offsets->members) {
        my $bam_readcount_entry = $reader->get_entry($vcf_entry->{chrom},
            $vcf_entry->{position} + $offset);
        if (defined $bam_readcount_entry) {
            $allele_hash{$offset} = $bam_readcount_entry->to_string;
        }
        else {
            $allele_hash{$offset} = ".";
        }
    }
    return encode(to_json(\%allele_hash, {canonical => 1}));
}

sub add_readcount_to_vcf_entry {
    my ($readcount_entry, $vcf_entry, $sample_idx) = @_;

    return;
}

sub encode {
    my $line = shift;

    $line =~ s/([\t])/?/g;
    $line =~ s/([:])/;/g;
    return $line;
}

sub decode {
    my $line = shift;

    $line =~ s/([?])/\t/g;
    $line =~ s/([;])/:/g;
    return $line;
}
1;

