package Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult::Index;

use strict;
use warnings;

use Genome;

#TODO - source directory
class Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::Index::Base',
    has => [
        config_file => {
            is => 'Text',
            is_calculated => 1,
            calculate => q{ "$output_dir/config.txt" },
            calculate_from => ["output_dir"],
        },
        files_index => {
            is => "Genome::Model::RnaSeq::DetectFusionsResult::FileIndex",
            is_transient => 1,
        }
    ],
    has_param => [
        ensembl_release => {
            is => 'Text',
            doc => 'the ensembl release number, used to determine the gtf file to download for use in deFuse'
        },
        config_params => {
            is => 'Text',
            doc => 'config params string, only params relevant to this index (and create_reference_dataset)'
        },
        version => {
            is => 'Text',
            doc => 'version of deFuse to run'
        }
    ],
    doc => 'create all the files needed for a successful run of deFuse'
};


sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/defuse-index/' . $self->id;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_) or return;

    unless ($self->reference_build->species_name eq "human"){
        die($self->error_message("tophat-fusion only supports human fusion detection"));
    }

    $self->_prepare_staging_directory;

    $self->files_index(Genome::Model::DetectFusionsResult::DeFuseResult::FileIndex->get_or_create(
        ensembl_release => $self->ensembl_release
    ));

    $self->generate_config_file;
    $self->_create_reference_dataset;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub _create_reference_dataset {
    my $self = shift;

    my $cmd = Genome::Model::RnaSeq::DetectFusionsResult->_path_for_command($self->version, "create_reference_dataset.pl");
    $cmd .= " -c " . $self->temp_staging_directory . "/config.txt";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->config_file],
        output_directories => [$self->temp_staging_directory . "/defuse"],
    );
}

sub _parse_params_into_hashref {
    my ($self, $params) = @_;
    #param syntax should be param_name=some value next_param_name=some other value
    #split into a key => value hash for use in building the deFuse config file
    my %params = $params =~ m/(\w+)=((?:[^\s=]+\s+)*[^\s=]+(?=\s+\w+=|$))/g;
    return \%params;
}

my %default_hash = (
    source_directory             => '[Where you unpacked the defuse code]',
    dataset_directory            => '[Where you intend to store the dataset]',
    gene_models                  => '[Gene models gtf file downloaded from ensembl or similar]',
    genome_fasta                 => '[Genome fasta]',
    repeats_filename             => '[Repeats table downloaded from ucsc]',
    est_fasta                    => '[EST fasta downloaded from ucsc]',
    est_alignments               => '[spliced EST alignments downloaded from ucsc]',
    unigene_fasta                => '[Unigene fasta downloaded from ncbi]',
    bowtie_bin                   => '[path of your bowtie binary version 0.12.5 or greater recommended]',
    bowtie_build_bin             => '[path of your bowtie-build binary version 0.12.5 or greater recommended]',
    blat_bin                     => '[path of your blat binary]',
    fatotwobit_bin               => '[path of your faToTwoBit binary from the blat suite]',
    r_bin                        => '[path of your R binary]',
    rscript_bin                  => '[path of your Rscript binary]',
    dataset_prefix               => '$(dataset_directory)/defuse',
    chromosome_prefix            => '$(dataset_prefix).dna.chromosomes',
    exons_fasta                  => '$(dataset_prefix).exons.fa',
    cds_fasta                    => '$(dataset_prefix).cds.fa',
    cdna_regions                 => '$(dataset_prefix).cdna.regions',
    cdna_fasta                   => '$(dataset_prefix).cdna.fa',
    reference_fasta              => '$(dataset_prefix).reference.fa',
    rrna_fasta                   => '$(dataset_prefix).rrna.fa',
    ig_gene_list                 => '$(dataset_prefix).ig.gene.list',
    repeats_regions              => '$(dataset_directory)/repeats.regions',
    est_split_fasta1             => '$(dataset_directory)/est.1.fa',
    est_split_fasta2             => '$(dataset_directory)/est.2.fa',
    est_split_fasta3             => '$(dataset_directory)/est.3.fa',
    est_split_fasta4             => '$(dataset_directory)/est.4.fa',
    est_split_fasta5             => '$(dataset_directory)/est.5.fa',
    est_split_fasta6             => '$(dataset_directory)/est.6.fa',
    est_split_fasta7             => '$(dataset_directory)/est.7.fa',
    est_split_fasta8             => '$(dataset_directory)/est.8.fa',
    est_split_fasta9             => '$(dataset_directory)/est.9.fa',
    prefilter1                   => '$(unigene_fasta)',
    scripts_directory            => '$(source_directory)/scripts',
    tools_directory              => '$(source_directory)/tools',
    data_directory               => '$(source_directory)/data',
    samtools_bin                 => '$(source_directory)/external/samtools-0.1.8/samtools',
    bowtie_threads               => '1',
    bowtie_quals                 => '--phred33-quals',
    max_insert_size              => '500',
    chromosomes                  => '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT',
    mt_chromosome                => 'MT',
    gene_sources                 => 'IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,processed_transcript,protein_coding',
    ig_gene_sources              => 'IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_pseudogene',
    rrna_gene_sources            => 'Mt_rRNA,rRNA,rRNA_pseudogene',
    num_blat_sequences           => '10000',
    dna_concordant_length        => '2000',
    discord_read_trim            => '50',
    clustering_precision         => '0.95',
    span_count_threshold         => '5',
    split_count_threshold        => '3',
    percent_identity_threshold   => '0.90',
    max_dist_pos                 => '600',
    num_dist_genes               => '500',
    split_min_anchor             => '4',
    max_concordant_ratio         => '0.1',
    splice_bias                  => '10',
    denovo_assembly              => 'no',
    positive_controls            => '$(data_directory)/controls.txt',
    probability_threshold        => '0.50',
    covariance_sampling_density  => '0.01',
    reads_per_job                => '1000000',
    regions_per_job              => '20',
    mailto                       => 'acoffman@genome.wustl.edu',
    remove_job_files             => 'yes',
    remove_job_temp_files        => 'yes',
    data_lane_regex_1            => '^(.+)_[12]_export\.txt.*$',
    data_end_regex_1             => '^.+_([12])_export\.txt.*$',
    data_compress_regex_1        => '^.+_[12]_export\.txt(.*)$',
    data_converter_1             => '$(scripts_directory)/fq_all2std.pl export2std',
    data_lane_regex_2            => '^(.+)_[12]_concat_qseq\.txt.*$',
    data_end_regex_2             => '^.+_([12])_concat_qseq\.txt.*$',
    data_compress_regex_2        => '^.+_[12]_concat_qseq\.txt(.*)$',
    data_converter_2             => '$(scripts_directory)/qseq2fastq.pl',
    data_lane_regex_3            => '^(.+)\.bam.*$',
    data_compress_regex_3        => '^.+\.bam(.*)$',
    data_end1_converter_3        => 'samtools view - | filter_sam_mate.pl 1 | sam_to_fastq.pl',
    data_end2_converter_3        => 'samtools view - | filter_sam_mate.pl 2 | sam_to_fastq.pl',
    data_lane_regex_4            => '^(.+).[12].fastq.*$',
    data_end_regex_4             => '^.+.([12]).fastq.*$',
    data_compress_regex_4        => '^.+.[12].fastq(.*)$',
    data_converter_4             => 'cat'
);

sub generate_config_file {
    my $self = shift;
    my $use_output_dir = shift;

    my $params_hash_ref = $self->_parse_params_into_hashref($self->config_params);

    #merge the user supplied params and the downloaded file locations into the default params
    my %config_hash = (
        %default_hash,
        %$params_hash_ref,
        dataset_directory => $use_output_dir ? $self->output_dir : $self->temp_staging_directory,
        log_directory => $self->temp_staging_directory,
        (map {$_ => $self->files_index->$_} qw{genome_fasta gene_models repeats_filename est_fasta est_alignments unigene_fasta}),
    );

    my $config_file = Genome::Sys->open_file_for_writing(($use_output_dir ? $use_output_dir : $self->temp_staging_directory) . "/config.txt");

    while(my ($key, $value) = each(%config_hash)){
        $config_file->say("$key = $value");
    }

    $config_file->close();

    $self->_reallocate_disk_allocation if $use_output_dir;
}


1;
