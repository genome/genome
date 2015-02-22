package Genome::Model::Tools::Vcf::AnnotateWithFpkm;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Writer;
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(firstidx);
use Memoize qw(memoize);
use Genome::File::Vcf::VepConsequenceParser;

my $FPKM_TAG = 'FPKM';
my $FPKM_HEADER = sprintf('<ID=%s,Number=A,Type=Float,Description="Framents Per Kilobase of exon per Million fragments">',$FPKM_TAG);

class Genome::Model::Tools::Vcf::AnnotateWithFpkm {
    is => 'Command::V2',
    has_input => [
        vcf_file => {
            is => 'File',
            doc => 'The vcf (gzipped or not) file that is to be annotated. Must be VEP annotated.',
        },
        fpkm_file => {
            is => 'File',
            doc => 'The file containing genes and fpkm values',
        },
        sample_name => {
            is => 'Text',
            doc => 'The sample name of the build containing the fpkm file',
        },
        output_file => {
            is_output => 1,
            doc => 'The output file, should end in ".vcf" or ".vcf.gz".  If ".vcf.gz" the output will be zipped using bgzip.',
        },
    ],
    has_transient_optional => [
        _transcripts_without_genes => {
            is => 'Hash',
        },
    ],
};
# TODO we probably want this information to be added to the sample field rather than info
# Travis brought this up, you'd have different information for tumor vs normal, we just happen to only care
# about tumor.
sub execute {
    my $self = shift;

    my $vcf_reader = Genome::File::Vcf::Reader->new($self->vcf_file);

    my $header = $vcf_reader->{header};
    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($header);
    $header->add_format_str($FPKM_HEADER);
    my $vcf_writer = Genome::File::Vcf::Writer->new($self->output_file, $header);
    my $sample_index = $header->index_for_sample_name($self->sample_name);

    while (my $entry = $vcf_reader->next) {
        $self->write_fpkm_for_sample($vep_parser, $entry, $sample_index);
        $vcf_writer->write($entry);
    }
    $self->_print_transcripts_without_genes;

    return 1;
}

sub write_fpkm_for_sample {
    my ($self, $vep_parser, $entry, $sample_index) = @_;

    my $gt = $entry->genotype_for_sample($sample_index);
    unless (defined($gt)) {
        $self->add_fpkm_to_entry($entry, '.', $sample_index);
        return;
    }

    my $gene_to_fpkm_map = $self->map_genes_to_fpkm;
    my @fpkms;
    for my $genotype_allele ($gt->get_alleles) {
        # Handle the reference allele which doesn't have an FPKM value
        if ($genotype_allele eq $entry->{reference_allele}) {
            push @fpkms, '.';
            next;
        }
        my ($transcript, @extra) = $vep_parser->transcripts($entry, $genotype_allele);

        if (not defined $transcript) {
            die $self->error_message("Vep returned no transcripts for genotype allele (%s) and vcf entry (%s)", $genotype_allele, $entry->to_string);
        } elsif (@extra) {
            die $self->error_message("Vep returned multiple transcripts, we want only one. Returned: %s", join(",", ( map{$_->{'feature'}} ($transcript, @extra) ) ) );
        }
        my $gene = $transcript->{'gene'};
        my $fpkm;
        if ($gene) {
            $fpkm = $gene_to_fpkm_map->{$gene};
        } else {
            $self->_add_transcript_without_gene($transcript);
            $fpkm = '.';
        }
        unless (defined $fpkm) {
            $self->status_message("Could not find a fpkm value for gene (%s) from fpkm file (%s).", $gene, $self->fpkm_file);
            $fpkm = '.';
        }
        push @fpkms, $fpkm;
    }
    $self->add_fpkm_to_entry($entry, join(',', @fpkms), $sample_index);
}

sub _add_transcript_without_gene {
    my ($self, $transcript_to_add) = @_;
    my $transcripts = $self->_transcripts_without_genes || {};
    $transcripts->{ $transcript_to_add->{consequence} }++;
    $self->_transcripts_without_genes($transcripts);
}

sub _print_transcripts_without_genes {
    my $self = shift;
    my $transcripts = $self->_transcripts_without_genes;
    $self->status_message("Transcripts by consequence type:");
    for my $consequence (keys %$transcripts) {
        $self->status_message("%s: %s", $consequence, $transcripts->{$consequence});
    }
}

# Return a hash mapping transcript names to fpkm values
sub map_genes_to_fpkm {
    my $self = shift;
    my $fh = Genome::Sys->open_file_for_reading($self->fpkm_file);

    $self->validate_header($fh);

    my %genes_to_fpkm;
    while (my $line = $fh->getline) {
        chomp $line;
        my @contents = split "\t", $line;
        my $gene_id = $contents[$self->gene_id_header_index];
        my $fpkm = $contents[$self->fpkm_header_index];
        my $fpkm_status = $contents[$self->fpkm_status_header_index];

        unless ($fpkm_status eq 'OK') {
            $genes_to_fpkm{$gene_id} = '.';
        }

        if (defined $genes_to_fpkm{$gene_id}) {
            next if $genes_to_fpkm{$gene_id} eq '.';
            $genes_to_fpkm{$gene_id} += $fpkm;
        }
        $genes_to_fpkm{$gene_id} = $fpkm;
    }

    return \%genes_to_fpkm;
}
memoize("map_genes_to_fpkm");

sub expected_fpkm_header {
    return qw(tracking_id class_code nearest_ref_id gene_id gene_short_name tss_id locus length coverage FPKM FPKM_conf_lo FPKM_conf_hi FPKM_status)
}

sub gene_id_header_index {
    my $self = shift;
    return firstidx{ $_ eq 'gene_id' } $self->expected_fpkm_header;
}

sub fpkm_header_index {
    my $self = shift;
    return firstidx{ $_ eq 'FPKM' } $self->expected_fpkm_header;
}

sub fpkm_status_header_index {
    my $self = shift;
    return firstidx{ $_ eq 'FPKM_status' } $self->expected_fpkm_header;
}

sub validate_header {
    my $self = shift;
    my $fh = shift;
    my $header = $fh->getline;
    chomp $header;
    unless ($header eq join("\t", $self->expected_fpkm_header) ) {
        die $self->error_message("Header found in file (%s) does not match our expectations.\nHeader:\n%s\nExpected:\n%s", $self->fpkm_file, $header, join("\t", $self->expected_fpkm_header) );
    }
}

sub add_fpkm_to_entry {
    my ($self, $entry, $fpkm, $sample_index) = @_;

    $entry->add_format_field("FPKM");
    my $existing_fpkm = $entry->sample_field($sample_index, "FPKM");
    if (defined($existing_fpkm)) {
        die $self->error_message("FPKM info field is already set on vcf entry at chrom/pos (%s %s), old value (%s) new value (%s)", $entry->{chrom}, $entry->{position}, $existing_fpkm, $fpkm);
    }
    $entry->set_sample_field($sample_index, "FPKM", $fpkm);
}

1;
