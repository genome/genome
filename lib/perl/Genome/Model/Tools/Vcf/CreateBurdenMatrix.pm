package Genome::Model::Tools::Vcf::CreateBurdenMatrix;

use Data::Dumper;
use Genome::File::Vcf::Reader;
use Genome::File::Vep::Reader;
use Genome::File::Vep::Writer;
use Genome;
use Sort::Naturally qw(ncmp nsort);
use List::Util qw(sum);

use strict;
use warnings;

class Genome::Model::Tools::Vcf::CreateBurdenMatrix {
    is => "Command::V2",
    has => [
        output_file => {
            is => 'Text',
            is_optional => 0,
            doc => "Output variant matrix format",
            is_input => 1,
            is_output => 1,
        },
        vcf_file => {
            is => 'Text',
            is_optional => 0,
            doc => "Merged Multisample Vcf containing mutations from all samples",
            is_input => 1,
        },
        sample_list_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Limit Samples in the Variant Matrix to Samples Within this File - Sample_Id should be the first column of a tab-delimited file, all other columns are ignored",
        },
        vep_annotation_file => {
            is => 'FilesystemPath',
            doc => 'Annotation file to determine which alleles meet requirements',
            is_input => 1,
        },
        burden_test_annotation_file => {
            is => 'Path',
            doc => 'Output annotation file of single or multiple annotation lines per variant',
            is_input => 1,
            is_output => 1,
        },
        consequence_types => {
            is => 'Text',
            doc => 'Which types of alterations to include in the matrix e.g.  ESSENTIAL_SPLICE_SITE,STOP_GAINED,STOP_LOST,NON_SYNONYMOUS_CODING',
            is_many => 1,
            default => ['ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST','NON_SYNONYMOUS_CODING'],
        },
    ],
    has_transient_optional => {
        _type_regex => {
            is => "Text",
        }
    },
};

sub _reformat_vep_entry {
    my $entry = shift;
    $entry->{gene} = $entry->{extra}{HGNC} if $entry->{extra}{HGNC};
    $entry->{uploaded_variation} = join("_", $entry->{gene}, $entry->{uploaded_variation});
    $entry->{uploaded_variation} =~ s/[^A-Za-z0-9]/_/g;
    return $entry;
}

# Given a list of sample names, return the indices of the samples we want to process as
# specified in the sample list file (or just all indices if no sample list file is given)
sub _sample_indices {
    my ($self, @sample_names) = @_;
    if ($self->sample_list_file) {
        my %limit;
        my $fh = Genome::Sys->open_file_for_reading($self->sample_list_file);
        my @names = map {my @x = split(/\t]/); $x[0]} <$fh>;
        chomp @names;
        @limit{@names} = ();
        return grep {exists $limit{$sample_names[$_]}} 0..$#sample_names;
    }
    return 0..$#sample_names;
}

# Given an array of Vep annotation entries, find those where the consequence matches the types
# of events we are interested in. Then, for each site (chr/pos/allele) that had an interesting
# event, grab ALL annotation for that site.
# The return value is a hash $h{gene}{alt} = [ entries... ]
sub _filter_annotation {
    my ($self, @entries) = @_;
    my @interesting_entries = grep {$_->{consequence} =~ $self->_type_regex} @entries;
    my %interesting_variants = map {join(" ", $_->{gene}, $_->{allele}) => 1} @interesting_entries;
    my %annotation;
    for my $entry (@entries) {
        my ($gene, $alt) = ($entry->{gene}, $entry->{allele});
        next unless exists $interesting_variants{"$gene $alt"};
        push(@{$annotation{$gene}{$alt}}, $entry);
    }
    return %annotation;
}

# Grab ALL annotation for every variant at the given chr/pos/@alts as long as at least one
# annotation entry matches one of the consequence types we are looking for.
sub _fetch_interesting_annotation {
    my ($self, $vep_in, $chr, $pos, @alts) = @_;
    my $loc = "$chr:$pos";
    while (my $p = $vep_in->peek) {
        last if ncmp($p->{location}, $loc) >= 0;
        $vep_in->next;
    }

    my @entries;
    while (my $p = $vep_in->peek) {
        last if ncmp($p->{location}, $loc) != 0;
        my $entry = $vep_in->next;
        next unless grep {$_ eq $entry->{allele}} @alts;
        push(@entries, _reformat_vep_entry($entry));
    }
    my %annotation = $self->_filter_annotation(@entries);
    return unless %annotation;

    return %annotation;
}

sub _variant_name {
    my ($gene, $vcf, @alts) = @_;

    # The original logic for computing the variant name is below.
    # I believe this corresponds to what Vep outputs.
    #my $id = '.';
    #$id = join("_", @{$vcf->{identifiers}}) if $vcf->{identifiers};
    #$id = join("_", $vcf->{chrom}, $vcf->{position}, join("_", $vcf->alleles)) if $id eq '.';
    #my $name = "${gene}_$id";
    #$name =~ s/[^A-Za-z0-9]/_/g;

    # We're now doing something a bit more descriptive.
    my @ids = grep {defined $_ && $_ ne '.'} @{$vcf->{identifiers}};
    my @fields = ($gene, @ids, $vcf->{chrom}, $vcf->{position}, $vcf->{reference_allele}, @alts);
    return join('_', @fields);
}

sub execute {
    my $self = shift;
    my $type_regex_str = join("|", $self->consequence_types);
    $self->_type_regex( qr/$type_regex_str/ );

    my $vcf_in = Genome::File::Vcf::Reader->new($self->vcf_file);
    my $vep_in = Genome::File::Vep::Reader->new($self->vep_annotation_file);
    my $vep_out = Genome::File::Vep::Writer->new($self->burden_test_annotation_file);
    my $mtx_out = Genome::Sys->open_file_for_writing($self->output_file);

    # Discover which samples we are interested in
    my @vcf_samples = $vcf_in->header->sample_names;
    my @sample_indices = $self->_sample_indices(@vcf_samples);
    $mtx_out->write(join("\t", "VariantId", @vcf_samples[@sample_indices]) . "\n");

    my $chrom_progress = ''; # For progress display
    # For each site
#my $count = 0;
    while (my $vcf_entry = $vcf_in->next) {
        print "Processing chromosome $vcf_entry->{chrom}\n" if $chrom_progress ne $vcf_entry->{chrom};
        $chrom_progress = $vcf_entry->{chrom};

        next if $vcf_entry->is_filtered; # We don't like sites that failed filters

        # Get total number of genotype calls (typically 2 per sample), and the associated
        # allele frequencies (counts, not relative frequency).
        my ($total, %allele_counts) = $vcf_entry->allelic_distribution;
        next unless $total; # no data for this site

        # If the reference is not present in any sample at this site, we will treat the major
        # allele as the reference. This means removing the major allele from the list of alleles
        # we are counting.
        my @alleles_by_desc_freq = sort { $allele_counts{$b} <=> $allele_counts{$a} } keys %allele_counts;
        shift @alleles_by_desc_freq if !exists $allele_counts{$vcf_entry->{reference_allele}};


        # Grab all interesting annotation for the alleles we found at the site.
        # The hash below is of the form $h{gene}{alt} = [annotation entries]
        my @all_alleles = $vcf_entry->alleles;
        my @allele_indices = keys %allele_counts;
        my @alleles_at_site = @all_alleles[@allele_indices];
        my %annotation = $self->_fetch_interesting_annotation(
            $vep_in, $vcf_entry->{chrom}, $vcf_entry->{position}, @alleles_at_site);
        next unless %annotation;

        # For each gene that had interesting annotation
        for my $gene (nsort keys %annotation) {
            $self->_process_gene($gene, \%annotation, \@sample_indices, $vcf_entry, $mtx_out, $vep_out);
        }
#exit if (++$count > 2000);
    }

    return 1;
}

sub _process_gene {
    my ($self, $gene, $annotation, $sample_indices, $vcf_entry, $mtx_out, $vep_out) = @_;

    my @alts_in_gene = sort keys %{$annotation->{$gene}};
    my @sample_alt_counts;
    my %alts_reported;
    my @all_alleles = $vcf_entry->alleles;

    # For each sample, count up the interesting alleles
    for my $si (@$sample_indices) {
        my ($total, %allele_idx_counts) = $vcf_entry->allelic_distribution($si);
        my @tmp = keys %allele_idx_counts;
        my %allele_counts; @allele_counts{@all_alleles[@tmp]} = @allele_idx_counts{@tmp};

        if ($total > 0) {
            my @interesting = grep {exists $allele_counts{$_} && $allele_counts{$_} > 0} @alts_in_gene;
            @alts_reported{@interesting} = 1;
            my $n = sum(map { $_ || 0 } @{allele_counts}{@interesting}) || 0;
            push(@sample_alt_counts, $n);
        } else {
            push(@sample_alt_counts, '.');
        }
    }

    # If we counted at least one thing (across all samples), we should tell the world about it.
    if (%alts_reported) {
        my $vname = _variant_name($gene, $vcf_entry, @alts_in_gene);
        $mtx_out->write(join("\t", $vname, @sample_alt_counts) . "\n");
        for my $alt (keys %alts_reported) {
            for my $entry (@{$annotation->{$gene}{$alt}}) {
                $entry->{uploaded_variation} = $vname;
                $vep_out->write($entry);
            }
        }
    }
}

1;
