package Genome::Model::Tools::Validation::LongIndelsPartTwo;

use warnings;
use strict;
use Genome;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Validation::LongIndelsPartTwo {
    is => 'Command',
    has_input => [
        output_dir => {
            is => 'String',
            doc => 'directory for output files from LongIndelsPartOne.pm; must contain a file "contigs.fa"',
        },
        tumor_val_model_copy_id => {
            is => 'Number',
            doc => 'model ID for copy of tumor validation model made in LongIndelsPartOne.pm',
        },
        normal_val_model_copy_id => {
            is => 'Number',
            doc => 'model ID for copy of normal validation model made in LongIndelsPartOne.pm',
        },
    ],
    has_optional_input => [
        contigs_file => {
            is => 'String',
            doc => 'contigs.fa file output from LongIndelsPartOne.pm',
            is_optional => 1,
        },
        tier_file_location => {
            is => 'String',
            doc => 'tiering file location to be used by gmt fast-tier fast-tier. Defaults to build36 location.',
            is_optional => 1,
            default => '/gscmnt/ams1100/info/model_data/2771411739/build102550711/annotation_data/tiering_bed_files_v3',
        },
    ],
    doc => 'Final steps in the validation of 3bp indels.',
};

sub help_detail {
    return <<EOS
    This tool performs the last steps (#'s 6-8) of the 3bp indel validation process outlined on this wiki page: https://gscweb.gsc.wustl.edu/wiki/Medical_Genomics/Nimblegen_Solid_Phase_Capture_Validation/Analysis#.3E3bp_Indels. It also then creates an annotate-able file using an adaptor 'gmt annotate adaptor indel-contig'. It also creates a bed file from this adapted list, and then uses fast-tier to tier the final calls. Lastly, the tool prints some details for possible manual review tickets as well, so be sure to SAVE THE STDOUT.
EOS
}

sub execute {

    my $self = shift;

    #parse input params
    my $output_dir = $self->output_dir;
    my $tumor_val_model_copy_id = $self->tumor_val_model_copy_id;
    my $normal_val_model_copy_id = $self->normal_val_model_copy_id;
    my $tier_file_location = $self->tier_file_location;

    #look for contigs.fa file
    my $contigs_file = $self->contigs_file;
    unless (defined $contigs_file) {
        $contigs_file = $output_dir . "/contigs.fa";
    }
    unless (-s $contigs_file) {
        $self->error_message("Contigs file $contigs_file does not exist with size > 0.");
        return;
    }

    #get builds and bam files
    my $tumor_val_model_copy = Genome::Model->get($tumor_val_model_copy_id) or die "Could not find tumor model with id $tumor_val_model_copy_id.\n";
    my $normal_val_model_copy = Genome::Model->get($normal_val_model_copy_id) or die "Could not find normal model with id $normal_val_model_copy_id.\n";
    my $tumor_build = $tumor_val_model_copy->last_succeeded_build or die "Could not find last succeeded build from tumor model $tumor_val_model_copy_id.\n";
    my $normal_build = $normal_val_model_copy->last_succeeded_build or die "Could not find last succeeded build from normal model $normal_val_model_copy_id.\n";
    my $tumor_bam = $tumor_build->whole_rmdup_bam_file or die "Cannot find tumor .bam.\n";
    my $normal_bam = $normal_build->whole_rmdup_bam_file or die "Cannot find normal .bam.\n";

    #get readcounts from the tumor sample
    my $tumor_rc_file = $output_dir . "/tumor_counts.mm80.csv";
    my $tumor_rc_cmd = Genome::Model::Tools::Validation::CountContigs->create(
        bam_file => $tumor_bam,
        contig_fasta_file => $contigs_file,
        maximum_mismatch_quality_sum => '80',
        output_file => $tumor_rc_file,
    );
    unless($tumor_rc_cmd->execute) {
        die "Failed to obtain tumor readcounts.\n";
    }

    #get readcounts from the normal sample
    my $normal_rc_file = $output_dir . "/normal_counts.mm80.csv";
    my $normal_rc_cmd = Genome::Model::Tools::Validation::CountContigs->create(
        bam_file => $normal_bam,
        contig_fasta_file => $contigs_file,
        maximum_mismatch_quality_sum => '80',
        output_file => $normal_rc_file,
    );
    unless ($normal_rc_cmd->execute) {
        die "Failed to obtain normal readcounts.\n";
    }

    #combine the readcounts to do somatic calling
    my $calls_file = $output_dir . "/combined_counts.csv";
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
    my $adaptor_cmd = Genome::Model::Tools::Annotate::Adaptor::IndelContig->create(
        contig_count_file => $somatic_calls_file,
        contig_fasta_file => $contigs_file,
        output_file => $adapted_somatic_calls,
    );
    unless ($adaptor_cmd->execute) {
        die "Failed to adapt the indels (using gmt annotate adaptor indel-contig).\n";
    }

    #convert annotation file to bed format for fast tiering
    my $somatic_bed_file = $adapted_somatic_calls . ".bed";
    my $convert_cmd = Genome::Model::Tools::Bed::Convert::Indel::AnnotationToBed->create(
        source => $adapted_somatic_calls,
        output => $somatic_bed_file,
    );
    unless ($convert_cmd->execute) {
        die "Failed to convert to .bed format.\n";
    }

    #fast tier bed format indel list
    my $tier1_calls = $somatic_bed_file . ".tier1";
    my $tiering_cmd = Genome::Model::Tools::FastTier::FastTier->create(
        variant_bed_file => $somatic_bed_file,
        tier_file_location => $tier_file_location,
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
    print "Tumor BAM with contigs: $tumor_bam\n";
    print "Normal BAM with contigs: $normal_bam\n";

    return 1;
}

1;
