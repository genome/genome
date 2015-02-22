package Genome::VariantReporting::Suite::FlankingRegions::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;
use IPC::Run qw(run);
use File::Basename qw(dirname);
use Genome::File::Vcf::Writer;
use Genome::File::Vcf::Reader;
use Memoize qw();


class Genome::VariantReporting::Suite::FlankingRegions::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        reference_fasta_lookup => {
            is => 'Text',
        },
        tumor_sample_name => {is => 'Text'},
    ],
    has_param => [
        flank_size => {is => 'Number', },
    ],
    has_transient => [
        reference_fasta => {is => 'Path'},
    ],
};

sub output_filename_base {
    return 'flanking_regions.vcf';
}

sub output_filename {
    my $self = shift;
    return $self->output_filename_base.'.gz';
}

sub _run {
    my $self = shift;

    $self->_annotate;
    $self->_zip;
    $self->status_message("Successfully ran");

    return;
}

sub _annotate {
    my $self = shift;
    $self->status_message("Reading from ".$self->input_vcf);
    my $in = Genome::File::Vcf::Reader->new($self->input_vcf);
    my $header = $in->header;
    my $alt_header_string = sprintf('<ID=FSAF,Number=A,Type=String,Description="%s bp flanking sequence to the left and right of variant">', $self->flank_size);
    my $ref_header_string = sprintf('<ID=FSRF,Number=A,Type=String,Description="%s bp flanking sequence to the left and right of reference allele">', $self->flank_size);
    $header->add_info_str($alt_header_string);
    $header->add_info_str($ref_header_string);
    print STDERR "Writing to ".$self->final_vcf_file."\n";
    my $out = Genome::File::Vcf::Writer->new($self->final_vcf_file, $header);

    while (my $variant = $in->next) {
        my $ft = $variant->sample_field($self->sample_index($header), "FT");
        unless ($variant->is_filtered or (defined $ft and ($ft ne '.' and $ft ne 'PASS'))) {
            my $left_sequence = $self->get_left_flanking_region($variant);
            my $right_sequence = $self->get_right_flanking_region($variant);
            $variant->info->{FSRF} = $left_sequence.$variant->{reference_allele}.$right_sequence;
            my @alt_sequences;
            for my $alt (@{$variant->{alternate_alleles}}) {
                push @alt_sequences, $left_sequence.$alt.$right_sequence;
            }
            $variant->info->{FSAF} = join(',', @alt_sequences);
        }
        $out->write($variant);
    }
    print STDERR "Done writing to ".$self->final_vcf_file."\n";
    $out->close;
}

sub sample_index {
    my $self = shift;
    my $header = shift;

    return $header->index_for_sample_name($self->tumor_sample_name);
}
Memoize::memoize("sample_index", LIST_CACHE => 'MERGE');

sub get_left_flanking_region {
    my $self = shift;
    my $variant = shift;
    my $chr = $variant->{chrom};
    my $start = $variant->{position} - $self->flank_size;
    if ($start < 0) {
        $start = 0;
    }
    my $end = $variant->{position} - 1;
    return $self->extract_sequence($chr, $start, $end);
}

sub get_right_flanking_region {
    my $self = shift;
    my $variant = shift;
    my $chr = $variant->{chrom};
    my $length = length($variant->{reference_allele});
    my $start = $variant->{position} + $length;
    if ($start < 0) {
        $start = 0;
    }
    my $end = $start + $self->flank_size - 1;
    return $self->extract_sequence($chr,$start,$end);
}

sub extract_sequence {
    my $self = shift;
    my ($chr, $start, $end) = @_;
    my $cmd =  sprintf("samtools faidx %s $chr:$start-$end | grep -v \">\"", $self->reference_fasta);
    my $sequence = `$cmd`;
    chomp $sequence;
    return $sequence;
}

sub final_vcf_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename_base);
}

sub _zip {
    my $self = shift;

    # deletes $self->final_vcf_file and creates $self->final_output_file
    run(['bgzip', $self->final_vcf_file]);
}

sub final_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename);
}
