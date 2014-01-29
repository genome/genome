package Genome::Model::Tools::Music::Survival;

use warnings;
use strict;
use Carp;
use Genome;
use IO::File;
use POSIX qw( WIFEXITED );

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Survival {
    is => 'Genome::Model::Tools::Music::Base',
    has_input => [
        bam_list => {
            is => 'Text',
            doc => "List of sample names to be included in the analysis. (See Description)",
        },
        maf_file => {
            is => 'Text',
            is_optional => 1,
            doc => "List of mutations in MAF format",
        },
        output_dir => {
            is => 'Text',
            is_output => 1,
            doc => "Directory where output files will be written",
        },
        genetic_data_type => {
            is => 'Text',
            is_optional => 1,
            default => "gene",
            doc => "Correlate clinical data to \"gene\" or \"variant\" level data",
        },
        numeric_clinical_data_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Table of samples (y) vs. numeric clinical data category (x)",
        },
        categorical_clinical_data_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Table of samples (y) vs. categorical clinical data category (x)",
        },
        glm_clinical_data_file => {
            is => 'Text',
            is_optional => 1,
            doc => "Clinical traits, mutational profiles, other mixed clinical data (See DESCRIPTION).",
        },
        phenotypes_to_include => {
            is => 'Text',
            is_optional => 1,
            doc => "Include only these genes and/or phenotypes in the anlaysis. (COMMA-DELIMITED)",
        },
        legend_placement => {
            is => 'Text',
            is_optional => 1,
            default => 'bottomleft',
            doc => "Choose one of 'bottomleft', 'topleft', 'topright', or 'bottomright'.",
        },
        skip_non_coding => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => "Skip non-coding mutations from the provided MAF file",
        },
        skip_silent => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => "Skip silent mutations from the provided MAF file",
        },
    ],
    doc => "Create survival plots and P-values for clinical and mutational phenotypes.",
};

sub help_synopsis {
    return <<HELP
 ... music survival \\
        --bam-list /path/myBamList.tsv \\
        --maf-file /path/myMAF.tsv \\
        --numeric-clinical-data-file /path/myNumericData.tsv \\
        --categorical-clinical-data-file /path/myClassData.tsv \\
        --output-dir /path/output_directory

 ... music survival \\
        --bam-list /path/myBamList.tsv \\
        --maf-file /path/myMAF.tsv \\
        --glm-clinical-data-file /path/myGLMClinicalData.tsv \\
        --output-dir /path/output_directory

 ... music survival \\
        --bam-list /path/myBamList.tsv \\
        --maf-file /path/myMAF.tsv \\
        --genetic-data-type 'gene' \\
        --glm-clinical-data-file /path/myGlmClinicalData.tsv \\
        --phenotypes-to-include 'Race,Gender,TP53' \\
        --output-dir /path/output_directory

HELP
}

sub help_detail {
    return <<HELP

This command performs survival analysis and plots survival curves for mutational data, as well as any clinical traits of interest as specified via the --phenotypes-to-include input parameter. The analyses performed include the Kaplan-Meier estimator followed by the Cox Proportional Hazards model. Outputs for each gene/clinical trait analyzed include survival curves, a hazard ratio (with confidence intervals), and P-values and FDRs describing the significance of the difference between survivors and non-survivors.

All clinical data files are searched for the required (case insensitive) "vital_status" and "days_to_last_followup" columns which are paired to phenotypes via sample IDs for the survival analysis. The first column of all clinical data files MUST contain the sample IDs, same as in other MuSiC tools. By default, analysis is performed on every gene present in the MAF. Optionally, the analysis may be limited to only specific genes by listing them (comma delimited) after the --phenotypes-to-include input parameter. Survival analysis may also be performed on other columns in the clinical data file by adding the column headers to the list of entries specified after the --phenotypes-to-include input parameter.

Here are some general guildelines for creating clinical data input files:

=over 4

=item * Headers are required.

=item * The first column of each clinical data file must contain sample IDs which match those in both the --bam-list and the MAF variant list (in the MAF, this is the Tumor_Sample_Barcode column, specifically).

=item * In at least one of the clinical data files input, columns with headers "vital_status" and "days_to_last_followup" (case insensitive) must exist. "vital_status" must be delineated by 1's and 0's, where 0 denotes 'living', and 1 denotes 'deceased'.

=back

Note that all input files must be tab-separated.

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
EOS
}

sub execute {

    # parse input arguments
    my $self = shift;
    my $bam_list = $self->bam_list;
    my $genetic_data_type = $self->genetic_data_type;
    my $legend_placement = $self->legend_placement;

    # handle phenotype inclusions
    my @phenotypes_to_include;
    my @clinical_phenotypes_to_include;
    my @mutated_genes_to_include;
    if ($self->phenotypes_to_include) { @phenotypes_to_include = split /,/,$self->phenotypes_to_include; }

    # check genetic data type
    unless ($genetic_data_type =~ /^gene|variant$/i) {
        $self->error_message("Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter.");
        return;
    }

    # load clinical data and analysis types
    my %clinical_data;
    if ($self->numeric_clinical_data_file) {
        $clinical_data{'numeric'} = $self->numeric_clinical_data_file;
    }
    if ($self->categorical_clinical_data_file) {
        $clinical_data{'categ'} = $self->categorical_clinical_data_file;
    }
    if ($self->glm_clinical_data_file) {
        $clinical_data{'glm'} = $self->glm_clinical_data_file;
    }

    # create array of all sample names possibly included from clinical data and MAF
    my @all_sample_names; # names of all the samples, no matter if they are mutated or not
    my $sampleFh = IO::File->new( $bam_list ) or die "Couldn't open $bam_list. $!\n";
    while( my $line = $sampleFh->getline )
    {
        next if ( $line =~ m/^#/ );
        chomp( $line );
        my ( $sample ) = split( /\t/, $line );
        push( @all_sample_names, $sample );
    }
    $sampleFh->close;

    # loop through clinical data files and assemble survival data hash (vital_status and days_to_last_followup required);
    my %survival_data;
    my $vital_status_flag = 0;
    my $days_to_last_follow_flag = 0;

    for my $clin_file (keys %clinical_data) {

        #check filehandle
        my $clin_fh = new IO::File $clinical_data{$clin_file},"r";
        unless ($clin_fh) {
            $self->error_message("Failed to open $clinical_data{$clin_file} for reading: $!");
            return;
        }

        #initiate variables to hold column info
        my %phenotypes_to_print;
        my $vital_status_col = 0;
        my $days_to_last_follow_col = 0;

        #parse header and record column locations for needed data
        my $header = $clin_fh->getline;
        my @header_fields = split /\t/,$header;
        for (my $i = 1; $i <= $#header_fields; $i++) { #sample ID should be in first column of file
            my $field = $header_fields[$i];
            if ($field =~ /vital_status|vitalstatus/i) { $vital_status_col = $i; $vital_status_flag++; }
            if ($field =~ /days_to_last_(follow_up|followup)|daystolastfollowup/i) { $days_to_last_follow_col = $i; $days_to_last_follow_flag++; }
            if (scalar grep { /^$field$/i } @phenotypes_to_include) { $phenotypes_to_print{$field} = $i; }
        }

        #read through clinical data file and store needed data in a hash
        while (my $line = $clin_fh->getline) {
            chomp $line;
            my @fields = split /\t/,$line;
            my $sample = $fields[0];
            unless (scalar grep { m/^$sample$/ } @all_sample_names) {
                $self->debug_message("Skipping sample $sample. (Sample is not in --bam-list).");
                next;
            }
            if ($vital_status_col) {
                my $vital_status;
                if ($fields[$vital_status_col] =~ /^(0|living)$/i) { $vital_status = 0; }
                elsif ($fields[$vital_status_col] =~ /^(1|deceased)$/i) { $vital_status = 1; }
                else { $vital_status = "NA"; }
                $survival_data{$sample}{'vital_status'} = $vital_status;
            }
            if ($days_to_last_follow_col) { $survival_data{$sample}{'days'} = $fields[$days_to_last_follow_col]; }
            for my $pheno (keys %phenotypes_to_print) { $survival_data{$sample}{$pheno} = $fields[$phenotypes_to_print{$pheno}]; }
        }
        $clin_fh->close;

        # record phenotypes included from clinical data
        push @clinical_phenotypes_to_include, keys %phenotypes_to_print;

    }

    # check for necessary header fields
    unless ($vital_status_flag) {
        $self->error_message('Clinical data does not seem to contain a column labeled "vital_status".');
        return;
    }
    unless ($days_to_last_follow_flag) {
        $self->error_message('Clnical data does not seem to contain a column labeled "days_to_last_followup".');
        return;
    }

    # create temporary files for R command
    my $survival_data_file = Genome::Sys->create_temp_file_path();
    my $mutation_matrix = Genome::Sys->create_temp_file_path();

    # print survival data (temp file)
    my $surv_fh = new IO::File $survival_data_file,"w" or die "Couldn't open survival data filehandle.";
    print $surv_fh join("\t","Sample","Days_To_Last_Followup","Vital_Status");
    if (@clinical_phenotypes_to_include) { print $surv_fh "\t" . join("\t",@clinical_phenotypes_to_include); }
    print $surv_fh "\n";
    for my $sample (keys %survival_data) {
        unless (exists $survival_data{$sample}{'days'}) { $survival_data{$sample}{'days'} = "NA"; }
        unless (exists $survival_data{$sample}{'vital_status'}) { $survival_data{$sample}{'vital_status'} = "NA"; }
        print $surv_fh join("\t",$sample,$survival_data{$sample}{'days'},$survival_data{$sample}{'vital_status'});
        for my $pheno (@clinical_phenotypes_to_include) {
            unless (exists $survival_data{$sample}{$pheno}) { $survival_data{$sample}{$pheno} = "NA"; }
            print $surv_fh "\t" . $survival_data{$sample}{$pheno};
        }
        print $surv_fh "\n";
    }
    $surv_fh->close;

    # find if any of the "phenotypes_to_include" are genes, and if so, limit the MAF mutation matrix to those genes
    my %clinical_pheno_to_include;
    @clinical_pheno_to_include{@clinical_phenotypes_to_include} = ();
    for my $item (@phenotypes_to_include) {
        push @mutated_genes_to_include,$item unless exists $clinical_pheno_to_include{$item};
    }
    my $mutated_genes_to_include = \@mutated_genes_to_include;

    # create mutation matrix file
    if( $genetic_data_type =~ /^gene$/i ) {
        $self->create_sample_gene_matrix_gene( $mutation_matrix, $mutated_genes_to_include, @all_sample_names );
    }
    elsif( $genetic_data_type =~ /^variant$/i ) {
        $self->create_sample_gene_matrix_variant( $mutation_matrix, $mutated_genes_to_include, @all_sample_names );
    }
    else {
        $self->error_message( "Please enter either \"gene\" or \"variant\" for the --genetic-data-type parameter." );
        return;
    }

    # check and prepare output directory
    my $output_dir = $self->output_dir . "/";
    unless (-e $output_dir) {
        $self->debug_message("Creating output directory: $output_dir...");
        unless(mkdir $output_dir) {
            $self->error_message("Failed to create output directory: $!");
            return;
        }
    }

    # set up R command
    my $R_cmd = "R --slave --args < " . __FILE__ . ".R " . join( " ", $survival_data_file, $mutation_matrix, $legend_placement, $output_dir );
    print "R_cmd:\n$R_cmd\n";

    #run R command
    WIFEXITED( system $R_cmd ) or croak "Couldn't run: $R_cmd ($?)";

    return(1);
}

sub create_sample_gene_matrix_gene {

    my ( $self, $mutation_matrix, $mutated_genes_to_include, @all_sample_names ) = @_;

    #create hash of mutations from the MAF file
    my ( %mutations, %all_genes );

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

        #check to see if this gene is on the list (if there is a list at all)
        if( defined($mutated_genes_to_include) and @{$mutated_genes_to_include} ) {
            next unless( scalar grep { m/^$gene$/ } @{$mutated_genes_to_include} );
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

    #sort @all_genes for consistency in header and loops
    my @all_genes = sort keys %all_genes;

    #write the input matrix for R code to a temp file
    my $matrix_fh = new IO::File $mutation_matrix,"w" or die "Failed to create matrix file $mutation_matrix!: $!";

    #print input matrix file header
    my $header = join( "\t", "Sample", @all_genes );
    $matrix_fh->print( "$header\n" );

    #print mutation relation input matrix
    for my $sample ( sort @all_sample_names ) {
        $matrix_fh->print( $sample );
        for my $gene ( @all_genes ) {
            if( exists $mutations{$sample}{$gene} ) {
                $matrix_fh->print( "\tMutated" );
            }
            else {
                $matrix_fh->print( "\tWildtype" );
            }
        }
        $matrix_fh->print( "\n" );
    }
}

sub create_sample_gene_matrix_variant {

    my ( $self, $mutation_matrix, $mutated_genes_to_include, @all_sample_names ) = @_;

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

        #check to see if this gene is on the list (if there is a list at all)
        if( defined @{$mutated_genes_to_include} ) {
            next unless (scalar grep { m/^$gene$/ } @{$mutated_genes_to_include});
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
    my $matrix_fh = new IO::File $mutation_matrix,"w" or die "Failed to create matrix file $mutation_matrix!: $!";

    #print input matrix file header
    my $header = join("\t","Sample",@variant_names);
    $matrix_fh->print("$header\n");

    #print mutation relation input matrix
    for my $sample (sort @all_sample_names) {
        $matrix_fh->print($sample);
        for my $variant (@variant_names) {
            if (exists $variants_hash{$sample}{$variant}) {
                $matrix_fh->print("\t$variants_hash{$sample}{$variant}");
            }
            else {
                $matrix_fh->print("\t0");
            }
        }
        $matrix_fh->print("\n");
    }
}

1;
