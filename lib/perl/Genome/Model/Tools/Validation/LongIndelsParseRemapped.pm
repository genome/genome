package Genome::Model::Tools::Validation::LongIndelsParseRemapped;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Validation::LongIndelsParseRemapped {
    is => "Command::V2",
    has => [
        contigs_file => {
            is => 'String',
            doc => 'Contigs fasta from first part of LongIndels validation process',
        },
        tumor_bam => {
            is => 'String',
            doc => 'Path to input tumor bam file',
            is_optional => 1,
        },
        normal_bam => {
            is => 'String',
            doc => 'Path to input normal bam file',
            is_optional => 1,
        },
        output_dir => {
            is_output => 1,
            is_input => 1,
            is => 'String',
            doc => 'Location to place output files',
        },
        tier_file_location => {
            is => 'String',
            doc => 'Path to tier file',
        },
        skip => {
            is => 'Boolean',
            default => 0,
        },
    ],
};

sub execute {
    my $self = shift;

    if ($self->skip) {
        $self->warning_message("skip signal received, skipping");
        return 1;
    }

    unless($self->tumor_bam and $self->normal_bam) {
        $self->error_message("Tumor and normal bams must both be provided");
        return;
    }

    my $output_dir = $self->output_dir;

    #get readcounts from the tumor sample
    my $tumor_rc_file = $output_dir . "/tumor_counts.mm80.csv";
    $self->debug_message("Starting count contigs for tumor");
    my $tumor_rc_cmd = Genome::Model::Tools::Validation::CountContigs->create(
        bam_file => $self->tumor_bam,
        contig_fasta_file => $self->contigs_file,
        maximum_mismatch_quality_sum => '80',
        output_file => $tumor_rc_file,
    );
    unless($tumor_rc_cmd->execute) {
        die "Failed to obtain tumor readcounts.\n";
    }

    #get readcounts from the normal sample
    my $normal_rc_file = $output_dir . "/normal_counts.mm80.csv";
    $self->debug_message("Starting count contigs for normal");
    my $normal_rc_cmd = Genome::Model::Tools::Validation::CountContigs->create(
        bam_file => $self->normal_bam,
        contig_fasta_file => $self->contigs_file,
        maximum_mismatch_quality_sum => '80',
        output_file => $normal_rc_file,
    );
    unless ($normal_rc_cmd->execute) {
        die "Failed to obtain normal readcounts.\n";
    }

    #combine the readcounts to do somatic calling
    my $calls_file = $output_dir . "/combined_counts.csv";
    $self->debug_message("Starting combine counts");
    my $calls_cmd = Genome::Model::Tools::Validation::CombineCounts->create(
        count_files => join(",",$normal_rc_file,$tumor_rc_file),
        file_labels => join(",","normal","tumor"),
        somatic_comparisons => "normal=>tumor",
        output_file => $calls_file,
    );
    unless ($calls_cmd->execute) {
        die "Failed to execute combine-counts.\n";
    }

    #grep to obtain somatic calls (and the header)
    my $somatic_calls_file = $calls_file . ".somatic";
    my $cmd = "grep \'contig\\\|Somatic\' $calls_file > $somatic_calls_file";
    print `$cmd`;

    #create an annotate-able list of somatic calls, getting ref and var bases from the contigs
    #gmt annotate adaptor indel-contig --contig-count-file /gscmnt/sata204/info/medseq/prc1/validation/Indel_analysis/lcm3/long_indels/combined_counts_somatic.csv --contig-fasta-file /gscmnt/sata204/info/medseq/prc1/validation/Indel_analysis/lcm3/long_indels/contigs.fa --output-file /gscmnt/sata204/info/medseq/prc1/validation/Indel_analysis/lcm3/long_indels/combined_counts_somatic.csv.indel.contig
    my $adapted_somatic_calls = $somatic_calls_file . ".adapted";
    $self->debug_message("Running gmt annotate adaptor indel-contig");
    my $adaptor_cmd = Genome::Model::Tools::Annotate::Adaptor::IndelContig->create(
        contig_count_file => $somatic_calls_file,
        contig_fasta_file => $self->contigs_file,
        output_file => $adapted_somatic_calls,
    );
    unless ($adaptor_cmd->execute) {
        die "Failed to adapt the indels (using gmt annotate adaptor indel-contig).\n";
    }

    unless (-s $adapted_somatic_calls) {
        $self->warning_message("$adapted_somatic_calls is empty, which might be ok. Either way, there's no more work to be done!");
        return 1;
    }

    #convert annotation file to bed format for fast tiering
    my $somatic_bed_file = $adapted_somatic_calls . ".bed";
    $self->debug_message("Running gmt bed convert indel annotation-to-bed");
    my $convert_cmd = Genome::Model::Tools::Bed::Convert::Indel::AnnotationToBed->create(
        source => $adapted_somatic_calls,
        output => $somatic_bed_file,
    );
    unless ($convert_cmd->execute) {
        die "Failed to convert to .bed format.\n";
    }

    #fast tier bed format indel list
    my $tier1_calls = $somatic_bed_file . ".tier1";
    $self->debug_message("Running gmt fast-tier fast-tier");
    my $tiering_cmd = Genome::Model::Tools::FastTier::FastTier->create(
        variant_bed_file => $somatic_bed_file,
        tier_file_location => $self->tier_file_location,
        tier1_output => $tier1_calls,
    );
    unless ($tiering_cmd->execute) {
        die "Failed to fast-tier.\n";
    }

    #print details for a manual review ticket of tier 1 events
    my $wc = `wc -l $tier1_calls`;
    ($wc) = split /\s+/,$wc;
    print "\n\n\n############### IMPORTANT NOTES ##########################################\n\n";
    print "You may want to include these details in any manual review ticket submitted:\n\n";
    print "3bp Indels: $tier1_calls (" . $wc . " sites)\n";
    print "Tumor BAM with contigs: $self->tumor_bam\n";
    print "Normal BAM with contigs: $self->normal_bam\n";

    return 1;
}

1;

