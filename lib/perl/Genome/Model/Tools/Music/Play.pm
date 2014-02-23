package Genome::Model::Tools::Music::Play;

use strict;
use warnings;

use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Play {
    is => 'Command::V2',
    has_input => [
        bam_list => {
            is => 'Text',
            doc => 'Tab delimited list of BAM files [sample_name normal_bam tumor_bam]'
        },
        roi_file => {
            is => 'Text',
            doc => 'Tab delimited list of ROIs [chr start stop gene_name]'
        },
        reference_sequence => {
            is => 'Text',
            doc => 'Path to reference sequence in FASTA format'
        },
        output_dir => {
            is => 'Text', is_output => 1,
            doc => 'Directory where output files and subdirectories will be written',
        },
        maf_file => {
            is => 'Text',
            doc => 'List of mutations using TCGA MAF specifications v2.3'
        },
        pathway_file => {
            is => 'Text',
            doc => 'Tab-delimited file of pathway information',
        },
    ],
    has_optional_input => [
        numeric_clinical_data_file => {
            is => 'Text',
            doc => 'Table of samples (y) vs. numeric clinical data category (x)',
        },
        categorical_clinical_data_file => {
            is => 'Text',
            doc => 'Table of samples (y) vs. categorical clinical data category (x)',
        },
        numerical_data_test_method => {
            is => 'Text', default => 'cor',
            doc => "Either 'cor' for Pearson Correlation or 'wilcox' for the Wilcoxon Rank-Sum Test for numerical clinical data.",
        },
        glm_model_file => {
            is => 'Text',
            doc => 'File outlining the type of model, response variable, covariants, etc. for the GLM analysis. (See DESCRIPTION).',
        },
        glm_clinical_data_file => {
            is => 'Text',
            doc => 'Clinical traits, mutational profiles, other mixed clinical data (See DESCRIPTION).',
        },
        use_maf_in_glm => {
            is => 'Boolean', default => 0,
            doc => 'Set this flag to use the variant matrix created from the MAF file as variant input to GLM analysis.',
        },
        omimaa_dir => {
            is => 'Path',
            doc => 'omim amino acid mutation database folder',
            default => Genome::Sys->dbpath('omim', 'latest'),
        },
        cosmic_dir => {
            is => 'Path',
            doc => 'cosmic amino acid mutation database folder',
            default => Genome::Sys->dbpath('cosmic', 'latest'),
        },
        verbose => {
            is => 'Boolean', default => 1,
            doc => 'turn on to display larger working output',
        },
        clinical_correlation_matrix_file => {
            is => 'Text',
            doc => 'Optionally store the sample-vs-gene matrix used internally during calculations.',
        },
        mutation_matrix_file => {
            is => 'Text',
            doc => 'Optionally store the sample-vs-gene matrix used during calculations.',
        },
        permutations => {
            is => 'Number',
            doc => 'Number of permutations used to determine P-values',
        },
        normal_min_depth => {
            is => 'Integer',
            doc => 'The minimum read depth to consider a Normal BAM base as covered',
        },
        tumor_min_depth => {
            is => 'Integer',
            doc => 'The minimum read depth to consider a Tumor BAM base as covered',
        },
        min_mapq => {
            is => 'Integer',
            doc => 'The minimum mapping quality of reads to consider towards read depth counts',
        },
        show_skipped => {
            is => 'Boolean', default => 0,
            doc => 'Report each skipped mutation, not just how many',
        },
        genes_to_ignore => {
            is => 'Text',
            doc => 'Comma-delimited list of genes to ignore for background mutation rates',
        },
        bmr => {
            is => 'Number',
            doc => 'Background mutation rate in the targeted regions',
        },
        max_proximity => {
            is => 'Text',
            doc => 'Maximum AA distance between 2 mutations',
        },
        bmr_modifier_file => {
            is => 'Text',
            doc => 'Tab delimited list of values per gene that modify BMR before testing [gene_name bmr_modifier]',
        },
        downsample_large_genes => {
            is => 'Boolean',
            default => 0,
            doc => "Downscale #bps in large genes, and the #muts proportionally"
        },
        skip_low_mr_genes => {
            is => 'Boolean', default => 1,
            doc => "Skip testing genes with MRs lower than the background MR"
        },
        max_fdr => {
            is => 'Number', default => 0.20,
            doc => 'The maximum allowed false discovery rate for a gene to be considered an SMG',
        },
        genetic_data_type => {
            is => 'Text',
            doc => 'Data in matrix file must be either "gene" or "variant" type data',
        },
        wu_annotation_headers => {
            is => 'Boolean',
            doc => 'Use this to default to wustl annotation format headers',
        },
        bmr_groups => {
            is => 'Integer', default => 1,
            doc => 'Number of clusters of samples with comparable BMRs',
        },
        separate_truncations => {
            is => 'Boolean', default => 0,
            doc => 'Group truncational mutations as a separate category',
        },
        merge_concurrent_muts => {
            is => 'Boolean', default => 0,
            doc => 'Multiple mutations of a gene in the same sample are treated as 1',
        },
        skip_non_coding => {
            is => 'Boolean', default => 1,
            doc => 'Skip non-coding mutations from the provided MAF file',
        },
        skip_silent => {
            is => 'Boolean', default => 1,
            doc => 'Skip silent mutations from the provided MAF file',
        },
        min_mut_genes_per_path => {
            is => 'Integer', default => 1,
            doc => 'Pathways with fewer mutated genes than this will be ignored',
        },
        processors => {
            is => 'Integer', default => 1,
            doc => "Number of processors to use in SMG (requires 'foreach' and 'doMC' R packages)",
        },
        aa_range => {
            is => 'Integer', default => 2,
            doc => "Set how close a 'near' match is when searching for amino acid near hits",
        },
        nuc_range => {
            is => 'Integer', default => 5,
            doc => "Set how close a 'near' match is when searching for nucleotide position near hits",
        },
        reference_build => {
            is => 'Text', default => 'Build37',
            doc => 'Put either "Build36" or "Build37"',
            example_values => ['Build36', 'Build37'],
        },
        show_known_hits => {
            is => 'Boolean', default => 1,
            doc => "When a finding is novel, show known AA in that gene",
        },
    ],
    has_calculated_optional => [
        gene_covg_dir => {
            calculate_from => ['output_dir'],
            calculate => q{ $output_dir . '/gene_covgs'; },
        },
        gene_mr_file => {
            calculate_from => ['output_dir'],
            calculate => q{ $output_dir . '/gene_mrs'; },
        },
        gene_list => {
            is => 'Text',
            doc => 'List of genes to test in B<genome-music-mutation-relation>(1), typically SMGs. (Uses output from running B<genome-music-smg>(1).)',
            calculate_from => ['output_dir'],
            calculate => q{ $output_dir . '/smg'; },
        },
        input_clinical_correlation_matrix_file => {
            is => 'Text', is_optional => 1,
            doc => "Instead of calculating this from the MAF, input the sample-vs-gene matrix used internally during calculations.",
        },
        mutation_relation_file => {
            is => 'Text', is_optional => 1,
            doc => 'Results of mutation-relation tool',
            calculate_from => ['output_dir'],
            calculate => q{ $output_dir . '/mutation_relation.csv'; },
        },
    ],
    has_constant => [
        cmd_list_file => { #If a workflow version of this tool is written, these parameters might be more useful
            is => 'Text', default_value => undef, is_optional => 1,
        },
        cmd_prefix => {
            is => 'Text', default_value => undef, is_optional => 1,
        },
    ],
    doc => 'Run the full suite of MuSiC tools sequentially.',
};

sub help_synopsis {
    return <<EOS
This tool takes as parameters all the information required to run the individual tools. An example usage is:

... music play \\
    --bam-list input/bams_to_analyze.txt \\
    --numeric-clinical-data-file input/numeric_clinical_data.csv \\
    --maf-file input/myMAF.tsv \\
    --output-dir play_output_dir \\
    --pathway-file input/pathway_db \\
    --reference-sequence input/refseq/all_sequences.fa \\
    --roi-file input/all_coding_regions.bed \\
    --genetic-data-type gene
EOS
}

sub help_detail {
    return <<EOS
This command can be used to run all of the MuSiC analysis tools on a set of data. Please see the individual tools for further description of the parameters.
EOS
}

sub _doc_credits {
    return "Please see the credits for B<genome-music>(1).";
}

sub _doc_authors {
    return " Thomas B. Mooney, M.S.";
}

sub _doc_see_also {
    return <<EOS
B<genome-music>(1),
B<genome-music-path-scan>(1),
B<genome-music-smg>(1),
B<genome-music-clinical-correlation>(1),
B<genome-music-mutation-relation>(1),
B<genome-music-cosmic-omim>(1),
B<genome-music-proximity>(1),
B<genome-music-pfam>(1)
EOS
}

sub execute {
    my $self = shift;

    my @no_dependencies = ('Proximity', 'ClinicalCorrelation', 'CosmicOmim', 'Pfam');
    my @bmr = ('Bmr::CalcCovg', 'Bmr::CalcBmr');
    my @depend_on_bmr = ('PathScan', 'Smg');
    my @depend_on_smg = ('MutationRelation');
    my @depends_on_all_others = ('CreateVisualizations'); #TODO: if this requires new params, insert them above
    for my $command_name (@no_dependencies, @bmr, @depend_on_bmr, @depend_on_smg, @depends_on_all_others) {
        my $command = $self->_create_command($command_name)
            or return;

        $self->_run_command($command)
            or return;
    }

    return 1;
}

sub _create_command {
    my $self = shift;
    my $command_name = shift;

    my $command_module = join('::', 'Genome::Model::Tools::Music', $command_name);
    my $command_meta = $command_module->__meta__;

    my %params;
    for my $property ($command_meta->_legacy_properties()) {
        next unless exists $property->{is_input} and $property->{is_input};
        my $property_name = $property->property_name;
        if($property_name eq 'output_file') {
            $params{$property_name} = $self->output_dir . '/' . $command_module->command_name_brief;
        } elsif(!$property->is_optional or defined $self->$property_name) {
            $params{$property_name} = $self->$property_name;
        }
    }

    my $command = $command_module->create(%params);
    unless($command) {
        $self->error_message('Failed to create command for ' . $command_name);
        return;
    }

    return $command;
}

sub _run_command {
    my $self = shift;
    my $command = shift;

    my $command_name = $command->command_name;
    $self->debug_message('Running ' . $command_name . '...');
    my $rv = eval { $command->execute() };
    if($@) {
        my $error = $@;
        $self->error_message('Error running ' . $command_name . ': ' . $error);
        return;
    } elsif(not $rv) {
        $self->error_message('Command ' . $command_name . ' did not return a true value.');
        return;
    } else {
        $self->debug_message('Completed ' . $command_name . '.');
        return 1;
    }
}

1;
