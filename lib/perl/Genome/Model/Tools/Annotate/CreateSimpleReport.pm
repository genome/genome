package Genome::Model::Tools::Annotate::CreateSimpleReport;

use warnings;
use strict;

use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::VepConsequenceParser;
use Genome::File::Vcf::Genotype;
use Genome::File::Tsv;

# gmt annotate create-simple-report
# --input-vcf-file=/gscuser/aregier/scratch/RT/AT-365/snvs.detailed.vep_segdup_hgvs.vcf.gz
# --sample-name=TEST-patient1-somval_tumor1
# --output-vcf-file=~/tmp/simple_report.tsv

class Genome::Model::Tools::Annotate::CreateSimpleReport {
    is => 'Command::V2',
    has_input => [
        input_vcf_file => {
            is => 'Text',
        },
        sample_name => {
            is => 'Text',
        },
        output_vcf_file => {
            is => 'Text',
        }
    ],
    has_transient_optional => [
        _sample_index =>  {
            is => 'Number',
            is_mutable => 0,
            calculate_from => ['input_vcf_file', 'sample_name'],
            calculate => q{
                my $input_vcf_fh = new Genome::File::Vcf::Reader($input_vcf_file);
                return $input_vcf_fh->header->index_for_sample_name($sample_name);
            },
        }
    ],
};

sub execute {
    my $self = shift;

    my $vcf_file = $self->input_vcf_file;
    my $sample_name = $self->sample_name;
    my $sample_index = $self->_sample_index;

    $self->validate_inputs;

    my $output_vcf_fh = Genome::Sys->open_file_for_appending($self->output_vcf_file);
    $self->print_header($output_vcf_fh);

    my $vcf = new Genome::File::Vcf::Reader($vcf_file);
    while (my $entry = $vcf->next) {
        next unless $entry->sample_field($sample_index, 'FT') eq 'PASS';

        my $variant_allele = $self->variant_allele($entry);
        next if $variant_allele eq $entry->{reference_allele};

        $self->print_basic_fields($output_vcf_fh, $entry, $variant_allele);
        print $output_vcf_fh "\t";
        $self->print_vep_fields($output_vcf_fh, $vcf->header, $entry, $variant_allele);
        print $output_vcf_fh "\n";
    }

    $output_vcf_fh->close;
    $self->validate_outputs;

    return 1;
}

#TODO this acutally better belongs in one of the Genome::File::Vcf classes
sub variant_allele {
    my $self = shift;
    my $entry = shift;

    my $genotype = $entry->genotype_for_sample($self->_sample_index);
    my @entry_allele_amino_acids = $entry->alleles;
    my @genotype_allele_pointers = $genotype->get_alleles;
    my @genotype_allele_amino_acids = map {$entry_allele_amino_acids[$_]} @genotype_allele_pointers;
    return $genotype_allele_amino_acids[1];
}

sub transcript {
    my $self = shift;
    my $vcf_header = shift;
    my $entry = shift;
    my $variant_allele = shift;

    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($vcf_header);
    my @transcripts = $vep_parser->transcripts($entry, $variant_allele);
    my $transcript;
    if (scalar(@transcripts) == 1) {
        $transcript = $transcripts[0];
    }
    elsif (scalar(@transcripts) > 1) {
        my @canonical_transcripts = $vep_parser->canonical_transcripts($entry, $variant_allele);
        if (scalar(@canonical_transcripts) == 1) {
            $transcript = $canonical_transcripts[0];
        }
        else {
            #TODO we need to handle this case
        }
    }
    else {
        #TODO we need to handle this case
    }

    return $transcript;
}

sub print_header {
    my $self = shift;
    my $output_vcf_fh = shift;

    my @headers = qw/
        chromosome_name
        start
        stop
        reference
        variant
        transcript_name
        trv_type
        amino_acid_change
        default_gene_name
        ensembl_gene_id
    /;

    print $output_vcf_fh join("\t", @headers) . "\n";
}

sub print_basic_fields {
    my $self = shift;
    my $output_vcf_fh = shift;
    my $entry = shift;
    my $variant_allele = shift;

    print $output_vcf_fh join("\t",
        $entry->{chrom},
        $entry->{position},
        $entry->{position} + length($variant_allele) - 1,
        $entry->{reference_allele},
        $variant_allele
    );
}

sub print_vep_fields {
    my $self = shift;
    my $output_vcf_fh = shift;
    my $vcf_header = shift;
    my $entry = shift;
    my $variant_allele = shift;

    my $transcript = $self->transcript($vcf_header, $entry, $variant_allele);
    unless (defined($transcript)) {
#TODO we need to handle this case
        print $output_vcf_fh "\tMore than one transcript found";
    }
    else {
        print $output_vcf_fh "\t" . join("\t",
            $transcript->{'feature'},
            $transcript->{'consequence'},
            $transcript->{'hgvsp'} || '',
            $transcript->{'symbol'} ,
            $transcript->{'gene'},
        );
    }
}

sub validate_inputs {
    my $self = shift;

    Genome::Sys->validate_file_for_reading($self->input_vcf_file);
    Genome::Sys->validate_file_for_writing($self->output_vcf_file);

    unless (defined($self->_sample_index)) {
        die $self->error_message(
            "Sample %s doesn't exist in input vcf %s",
            $self->sample_name,
            $self->input_vcf_file
        );
    }

    return 1;
}

sub validate_outputs {
    my $self = shift;

    unless (-s $self->output_vcf_file) {
        die $self->error_message("Output file %s doesn't exist or has non-zero size", $self->output_vcf_file);
    }

    return 1;
}

1;
