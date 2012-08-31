#!/usr/bin/env genome-perl
use strict;
use warnings;

use above 'Genome';
use Bio::DB::Sam;

#AML103
#my $rna_seq_model_id = '2880023893';
#my $som_var_model_id = '2879908905';

#HG1
my $rna_seq_model_id = '2876830807';
my $som_var_model_id = '2875816457';



my $annotation_output_version = 1;

my $som_var_model = Genome::Model->get($som_var_model_id);
my $som_var_build = $som_var_model->last_succeeded_build;

my $rna_seq_model = Genome::Model->get($rna_seq_model_id);
my $rna_seq_build = $rna_seq_model->last_succeeded_build;
unless ($rna_seq_build) {
    die('Failed to find RNA-seq build for model: '. $rna_seq_model_id);
}

# TODO: remove hard coded build
#my $rna_seq_build = Genome::Model::Build->get('114919792');

# Get the RNA-seq(probably TopHat) BAM file
my $alignment_result = $rna_seq_build->alignment_result;
my $rna_seq_bam_path = $alignment_result->bam_file;
print STDERR 'RNA-seq BAM file: '. $rna_seq_bam_path ."\n";

# Get the reference genome for the RNA-seq model
my $reference_build = $rna_seq_model->reference_sequence_build;
my $reference_path = $reference_build->full_consensus_path('fa');
print STDERR 'Reference Genome: '. $reference_path ."\n";

my $bam = Bio::DB::Bam->open($rna_seq_bam_path);
my $header = $bam->header;
my $index = Bio::DB::Bam->index($rna_seq_bam_path);
my $fai = Bio::DB::Sam::Fai->load($reference_path);

my $annotated_file = $som_var_build->data_set_path("effects/snvs.hq.tier1",$annotation_output_version,'annotated');

print STDERR 'Annotation File: '. $annotated_file ."\n";

# TODO: Figure out the remaining 7 column headers
my @headers = qw/
                    chr
                    start
                    stop
                    reference
                    variant
                    variation_type
                    gene
                    transcript
                    species
                    transcript_source
                    transcript_version
                    strand
                    transcript_status
                    type
                    unk_aa
                    unk_1
                    unk_2
                    unk_3
                    unk_4
                    unk_5
                    unk_6
                /;

my $reader = Genome::Utility::IO::SeparatedValueReader->create(
    headers => \@headers,
    input => $annotated_file,
    separator => "\t",
);

my $callback = sub {
    my ($tid,$pos,$pileups,$callback_data) = @_;
    my $data = $callback_data->[0];
    my $read_counts = $callback_data->[1];
    if ( ($pos == ($data->{start} - 1) ) ) {
        #print STDERR 'PILEUP:'. $data->{chr} ."\t". $tid ."\t". $pos ."\t". $data->{start} ."\t". $data->{stop}."\n";
        my $ref_base = $fai->fetch($data->{chr} .':'. $data->{start} .'-'. $data->{stop});
        unless ($data->{reference} eq $ref_base) {
            die('Reference base '. $ref_base .' does not match expected '. $data->{reference} .' at postion '. $pos .' for chr '. $data->{chr});
        }
        for my $pileup ( @{$pileups} ) {
            my $alignment = $pileup->alignment;
            next if $pileup->indel or $pileup->is_refskip;      # don't deal with these ;-)

            my $qbase  = substr($alignment->qseq,$pileup->qpos,1);
            next if $qbase =~ /[nN]/;
            $read_counts->{$qbase}++;
        }
    }
};

my @output_headers = (@headers,'A_count','T_count','G_count','C_count');

my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
    #output => 'tmp.tsv',
    headers => \@output_headers,
    separator => "\t",
    print_headers => 1,
);
while (my $data = $reader->next) {
    my $seq_id = $data->{chr} .':'. $data->{start} .'-'. $data->{stop};
    #print STDERR 'SeqId: '. $seq_id ."\n";
    unless ($data->{variation_type} eq 'SNP') { die('Only handles SNPs!'); }
    my ($tid,$start,$end) = $header->parse_region($seq_id);
    #print STDERR 'TargetId: '. $tid ."\n";
    #print STDERR 'Start: '. $start ."\n";
    #print STDERR 'End: '. $end ."\n";
    my %read_counts;
    $index->pileup($bam,$tid,$start,$end,$callback,[$data,\%read_counts]);
    $data->{A_count} = $read_counts{A} || 0;
    $data->{T_count} = $read_counts{T} || 0;
    $data->{C_count} = $read_counts{C} || 0;
    $data->{G_count} = $read_counts{G} || 0;
    $writer->write_one($data);
}


exit;
