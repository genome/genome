package Genome::Model::Tools::Music::ClinicalCorrelation;

use warnings;
use strict;
use Carp;
use Genome;
use IO::File;
use POSIX qw( WIFEXITED );

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::ClinicalCorrelation {
    is => 'Genome::Model::Tools::Music::Base',
    has_input => [
        bam_list => {
            is => 'Text',
            doc => "Tab delimited list of BAM files [sample_name, normal_bam, tumor_bam] (See Description)",
        },
        maf_file => {
            is => 'Text', is_optional => 1,
            doc => "List of mutations using TCGA MAF specification v2.3",
        },
        output_file => {
            is_output => 1, is => 'Text',
            doc => "Results of clinical-correlation tool. Will have suffix added for data type",
        },
        clinical_correlation_matrix_file => {
            is => 'Text', is_optional => 1,
            doc => "Specify a file to store the sample-vs-gene matrix created during calculations",
        },
        input_clinical_correlation_matrix_file => {
            is => 'Text', is_optional => 1,
            doc => "Instead of creating this from the MAF, input the sample-vs-gene matrix for calculations",
        },
        genetic_data_type => {
            is => 'Text', is_optional => 1, default => "gene",
            doc => "Correlate clinical data to \"gene\" or \"variant\" level data",
        },
        numeric_clinical_data_file => {
            is => 'Text', is_optional => 1,
            doc => "Table of samples (y) vs. numeric clinical data category (x)",
        },
        numerical_data_test_method => {
            is => 'Text', is_optional => 1, default => 'wilcox',
            doc => "Either 'cor' for Pearson Correlation or 'wilcox' for the Wilcoxon Rank-Sum Test for numerical clinical data",
        },
        categorical_clinical_data_file => {
            is => 'Text', is_optional => 1,
            doc => "Table of samples (y) vs. categorical clinical data category (x)",
        },
        glm_model_file => {
            is => 'Text', is_optional => 1,
            doc => "File outlining the type of model, response variable, covariants, etc. for the GLM analysis. (See DESCRIPTION)",
        },
        glm_clinical_data_file => {
            is => 'Text', is_optional => 1,
            doc => "Clinical traits, mutational profiles, other mixed clinical data (See DESCRIPTION)",
        },
        use_maf_in_glm => {
            is => 'Boolean', is_optional => 1, default => 0,
            doc => "Create a variant matrix from the MAF file as variant input to GLM analysis.",
        },
        skip_non_coding => {
            is => 'Boolean', is_optional => 1, default => 1,
            doc => "Skip non-coding mutations from the provided MAF file",
        },
        skip_silent => {
            is => 'Boolean', is_optional => 1, default => 1,
            doc => "Skip silent mutations from the provided MAF file",
        },
    ],
    doc => "Correlate phenotypic traits against mutated genes, or against individual variants",
};

sub help_synopsis {
    return <<HELP
 ... music clinical-correlation \\
        --bam-list /path/myBamList.tsv \\
        --maf-file /path/myMAF.tsv \\
        --numeric-clinical-data-file /path/myNumericData.tsv \\
        --genetic-data-type 'gene' \\
        --output-file /path/output_file

 ... music clinical-correlation \\
        --maf-file /path/myMAF.tsv \\
        --bam-list /path/myBamList.tsv \\
        --numeric-clinical-data-file /path/myNumericData.tsv \\
        --categorical-clinical-data-file /path/myClassData.tsv \\
        --genetic-data-type 'gene' \\
        --output-file /path/output_file

 ... music clinical-correlation \\
        --maf-file /path/myMAF.tsv \\
        --bam-list /path/myBamList.tsv \\
        --output-file /path/output_file \\
        --glm-model-file /path/model.tsv \\
        --glm-clinical-data-file /path/glm_clinical_data.tsv \\
        --use-maf-in-glm

HELP
}

sub help_detail {
    return <<HELP

This command relates clinical traits and mutational data. Either one can perform correlation analysis between mutations recorded in a MAF and the particular phenotypic traits recorded in clinical data files for the same samples, or one can run a generalized linear model (GLM) analysis on the same types of data.

The clinical data files for correlation must be separated between numeric and categoric data and must follow these
conventions:

=over 4

=item * Headers are required

=item * Each file must include at least 1 sample_id column and 1 attribute column, with the format being [sample_id  clinical_data_attribute_1 clinical_data_attribute_2 ...]

=item * The sample ID must match the sample ID listed in the MAF under "Tumor_Sample_Barcode" for relating the mutations of this sample.

=back

Note the importance of the headers: the header for each clinical_data_attribute will appear in the
output file to denote relationships with the mutation data from the MAF.

Internally, the input data is fed into an R script which calculates a P-value representing the
probability that the correlation seen between the mutations in each gene (or variant) and each
phenotype trait are random. Lower P-values indicate lower randomness, or likely true correlations.

The results are saved to the output filename given with a suffix appended; ".numeric.tsv" will be
appended for results derived from numeric clinical data, and ".categorical.tsv" will be appended for results
derived from categorical clinical data. Also, ".glm.tsv" will be appended to the output filename for GLM results.

The GLM analysis accepts a mixed numeric and categoric clinical data file, input using the parameter --glm-clinical-data-file. GLM clinical data must adhere to the formats described above for the correlation clinical data files. GLM also requires the user to input a --glm-model-file. This file requires specific headers and defines the analysis to be performed rather exactly. Here are the conventions required for this file:

=over 4

=item * Columns must be ordered as such:

=item [ analysis_type    clinical_data_trait_name    variant/gene_name   covariates  memo ]

=item * The 'analysis_type' column must contain either "Q", indicating a quantative trait, or "B", indicating a binary trait will be examined.

=item * The 'clinical_data_trait_name' is the name of a clinical data trait defined by being a header in the --glm-clinical-data-file.

=item * The 'variant/gene_name' can either be the name of one or more columns from the --glm-clinical-data-file, or the name of one or more mutated gene names from the MAF, separated by "|". If this column is left blank, or instead contains "NA", then each column from either the variant mutation matrix (--use-maf-in-glm) or alternatively the --glm-clinical-data-file is used consecutively as the variant column in independent analyses.

=item * 'covariates' are the names of one or more columns from the --glm-clinical-data-file, separated by "+".

=item * 'memo' is any note deemed useful to the user. It will be printed in the output data file for reference.

=back

GLM analysis may be performed using solely the data input into --glm-clinical-data-file, as described above, or alternatively, mutational data from the MAF may be included as variants in the GLM analysis, as also described above. Use the --use-maf-in-glm flag to include the mutation matrix derived from the maf as variant data.

Note that all input files for both correlation and GLM analysis must be tab-separated.

HELP
}

sub _additional_help_sections {
    return (
    "ARGUMENTS",
<<EOS

=over 4

=item --bam-list

=over 8

=item Provide a file containing sample names and normal/tumor BAM locations for each. Use the tab-
  delimited format [sample_name normal_bam tumor_bam] per line. This tool only needs sample_name,
  so all other columns can be skipped. The sample_name must be the same as the tumor sample names
  used in the MAF file (16th column, with the header Tumor_Sample_Barcode).

=back

=back

=cut

EOS
    );
}

sub _doc_authors {
    return <<EOS
 Nathan D. Dees, Ph.D.
 Qunyuan Zhang, Ph.D.
 William Schierding, M.S.
EOS
}

sub execute {

    # parse input arguments
    my $self = shift;
    my $bam_list = $self->bam_list;
    my $output_file = $self->output_file;
    my $genetic_data_type = $self->genetic_data_type;

    # check genetic data type
    unless( $genetic_data_type =~ /^gene|variant$/i ) {
        $self->error_message("Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter.");
        return;
    }

    # load clinical data and analysis types
    my %clinical_data;
    if( $self->numeric_clinical_data_file ) {
        $clinical_data{'numeric'} = $self->numeric_clinical_data_file;
    }
    if( $self->categorical_clinical_data_file ) {
        $clinical_data{'categ'} = $self->categorical_clinical_data_file;
    }
    if( $self->glm_clinical_data_file ) {
        $clinical_data{'glm'} = $self->glm_clinical_data_file;
    }
    my $glm_model = $self->glm_model_file;

    # declarations
    my @all_sample_names; # names of all the samples, no matter if it's mutated or not

    # parse out the sample names from the bam-list which should match the names in the MAF file
    my $sampleFh = IO::File->new( $bam_list ) or die "Couldn't open $bam_list. $!\n";
    while( my $line = $sampleFh->getline ) {
        next if ( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample ) = split( /\t/, $line );
        push( @all_sample_names, $sample );
    }
    $sampleFh->close;

    # loop through clinical data files
    for my $datatype ( keys %clinical_data ) {

        my $test_method;
        my $full_output_filename;

        if( $datatype =~ /numeric/i ) {
            $full_output_filename = $output_file . ".numeric.tsv";
            $test_method = $self->numerical_data_test_method;
        }

        if( $datatype =~ /categ/i ) {
            $full_output_filename = $output_file . ".categorical.tsv";
            $test_method = "fisher";
        }

        if( $datatype =~ /glm/i ) {
            $full_output_filename = $output_file . ".glm.tsv";
            $test_method = "glm";
        }

        #read through clinical data file to see which samples are represented and create input matrix for R
        my %samples;
        my $matrix_file;
        my $samples = \%samples;
        my $clin_fh = new IO::File $clinical_data{$datatype},"r";
        unless( $clin_fh ) {
            die "failed to open $clinical_data{$datatype} for reading: $!";
        }
        my $header = $clin_fh->getline;
        while( my $line = $clin_fh->getline ) {
            chomp $line;
            my ( $sample ) = split( /\t/, $line );
            $samples{$sample}++;
        }
        #create correlation matrix unless it's glm analysis without using a maf file
        unless(( $datatype =~ /glm/i && !$self->use_maf_in_glm ) || $self->input_clinical_correlation_matrix_file ) {

            if( $genetic_data_type =~ /^gene$/i ) {
                $matrix_file = $self->create_sample_gene_matrix_gene( $samples, $clinical_data{$datatype}, @all_sample_names );
            }
            elsif( $genetic_data_type =~ /^variant$/i ) {
                $matrix_file = $self->create_sample_gene_matrix_variant( $samples, $clinical_data{$datatype}, @all_sample_names );
            }
            else {
                $self->error_message( "Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter." );
                return;
            }
        }

        if( $self->input_clinical_correlation_matrix_file ) {
            $matrix_file = $self->input_clinical_correlation_matrix_file;
        }

        unless( defined $matrix_file ) { $matrix_file = "'*'"; }

        #set up R command
        my $R_cmd = "R --slave --args < " . __FILE__ . ".R $test_method ";

        if( $datatype =~ /glm/i ) {
            $R_cmd .= "$glm_model $clinical_data{$datatype} $matrix_file $full_output_filename";
        }
        else {
            $R_cmd .= "$clinical_data{$datatype} $matrix_file $full_output_filename";
        }

        #run R command
        print "R_cmd:\n$R_cmd\n";
        WIFEXITED( system $R_cmd ) or croak "Couldn't run: $R_cmd ($?)";
    }

    return( 1 );
}

sub create_sample_gene_matrix_gene {

    my ( $self, $samples, $clinical_data_file, @all_sample_names ) = @_;
    my $output_matrix = $self->clinical_correlation_matrix_file;

    #create a hash of mutations from the MAF file
    my ( %mutations, %all_genes, @all_genes );

    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $self->maf_file ) or die "Couldn't open MAF file!\n";
    while( my $line = $maf_fh->getline ) {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $mutation_class, $sample ) = @cols[0,8,15];

        #check that the mutation class is valid
        if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            $self->error_message( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }

        #check if sample exists in clinical data
        unless( defined $samples->{$sample} ) {
            warn "Sample Name: $sample from MAF file does not exist in Clinical Data File";
            next;
        }

        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(( $self->skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $self->skip_silent && $mutation_class =~ m/^Silent$/ )) {
            print "Skipping $mutation_class mutation in gene $gene.\n";
            next;
        }

        $all_genes{$gene}++;
        $mutations{$sample}{$gene}++;
    }
    $maf_fh->close;

    #sort @all_genes for consistency
    @all_genes = sort keys %all_genes;

    #write the input matrix for R code to a file
    my $matrix_fh;
    unless (defined $output_matrix) {
        $output_matrix = Genome::Sys->create_temp_file_path();
    }
    $matrix_fh = new IO::File $output_matrix,"w";
    unless ($matrix_fh) {
        die "Failed to create matrix file $output_matrix!: $!";
    }
    #print input matrix file header
    my $header = join("\t","Sample",@all_genes);
    $matrix_fh->print("$header\n");

    #print mutation relation input matrix
    for my $sample (sort @all_sample_names) {
        $matrix_fh->print($sample);
        for my $gene (@all_genes) {
            if (exists $mutations{$sample}{$gene}) {
                $matrix_fh->print("\t1");
            }
            else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }

    return $output_matrix;
}

sub create_sample_gene_matrix_variant {

    my ( $self, $samples, $clinical_data_file, @all_sample_names ) = @_;
    my $output_matrix = $self->clinical_correlation_matrix_file;

    #create hash of mutations from the MAF file
    my ( %variants_hash, %all_variants );

    #parse the MAF file and fill up the mutation status hashes
    my $maf_fh = IO::File->new( $self->maf_file ) or die "Couldn't open MAF file!\n";
    while( my $line = $maf_fh->getline ) {
        next if( $line =~ m/^(#|Hugo_Symbol)/ );
        chomp $line;
        my @cols = split( /\t/, $line );
        my ( $gene, $chr, $start, $stop, $mutation_class, $mutation_type, $ref, $var1, $var2, $sample ) = @cols[0,4..6,8..12,15];

        #check that the mutation class is valid
        if( $mutation_class !~ m/^(Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site|Translation_Start_Site|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Silent|Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|De_novo_Start_InFrame|De_novo_Start_OutOfFrame)$/ ) {
            $self->error_message( "Unrecognized Variant_Classification \"$mutation_class\" in MAF file for gene $gene\nPlease use TCGA MAF v2.3.\n" );
            return;
        }

        unless( exists $samples->{$sample} ) {
            warn "Sample Name: $sample from MAF file does not exist in Clinical Data File";
            next;
        }

        # If user wants, skip Silent mutations, or those in Introns, RNA, UTRs, Flanks, IGRs, or the ubiquitous Targeted_Region
        if(( $self->skip_non_coding && $mutation_class =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region)$/ ) ||
           ( $self->skip_silent && $mutation_class =~ m/^Silent$/ )) {
            print "Skipping $mutation_class mutation in gene $gene.\n";
            next;
        }

        my $var;
        my $variant_name;
        if( $ref eq $var1 ) {
            $var = $var2;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        }
        elsif( $ref eq $var2 ) {
            $var = $var1;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        }
        elsif( $ref ne $var1 && $ref ne $var2 ) {
            $var = $var1;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
            $var = $var2;
            $variant_name = $gene."_".$chr."_".$start."_".$stop."_".$ref."_".$var;
            $variants_hash{$sample}{$variant_name}++;
            $all_variants{$variant_name}++;
        }
    }
    $maf_fh->close;

    #sort variants for consistency
    my @variant_names = sort keys %all_variants;

    #write the input matrix for R code to a file
    my $matrix_fh;
    unless( defined $output_matrix ) {
        $output_matrix = Genome::Sys->create_temp_file_path();
    }
    $matrix_fh = new IO::File $output_matrix,"w";
    unless( $matrix_fh ) {
        die "Failed to create matrix file $output_matrix!: $!";
    }

    #print input matrix file header
    my $header = join( "\t", "Sample", @variant_names );
    $matrix_fh->print("$header\n");

    #print mutation relation input matrix
    for my $sample ( sort @all_sample_names ) {
        $matrix_fh->print( $sample );
        for my $variant ( @variant_names ) {
            if( exists $variants_hash{$sample}{$variant} ) {
                $matrix_fh->print("\t$variants_hash{$sample}{$variant}");
            }
            else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }

    return $output_matrix;
}

1;
