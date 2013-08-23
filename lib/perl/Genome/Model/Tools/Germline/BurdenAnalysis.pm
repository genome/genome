package Genome::Model::Tools::Germline::BurdenAnalysis;

use warnings;
use strict;
use Data::Dumper;
use Carp qw/confess/;
use File::Basename qw/dirname/;
use Genome;
use Workflow::Simple;
use Workflow;

my $base_R_commands = join("/", dirname(__FILE__), "gstat/burdentest/burdentest.R");

class Genome::Model::Tools::Germline::BurdenAnalysis {
    is => ['Genome::Command::Base'],
    has_input => [
        mutation_file => {
            is => 'Text',
            doc => "Mutation Matrix",
        },
        glm_clinical_data_file => {
            is => 'Text',
            doc => "Phenotype File",
        },
        VEP_annotation_file => {
            is => 'Text',
            doc => "List of mutations --VEP annotate then VEP parse",
        },
        project_name => {
            is => 'Text',
            doc => "The name of the project",
            default => 'VariantId',
        },
        output_directory => {
            is => 'Text',
            doc => "Results of the Burden Analysis",
            is_output => 1,
        },
        maf_cutoff => {
            is => 'Text',
            doc => "The cutoff to use to define which mutations are rare, 1 means no cutoff",
            default => '0.01'
        },
        permutations => {
            is => 'Text',
            doc => "The number of permutations to perform, typically from 100 to 10000, larger number gives more accurate p-value, but needs more time",
            default => '10000'
        },
        p_value_permutations => {
            is => 'Text',
            doc => "The number of permutations to perform on the p-values, typically from 10 to 100, larger number gives more accurate p-value, but needs more time",
            default => '100'
        },
        trv_types => {
            is => 'Text',
            doc => 'colon-delimited list of which trv types to use as significant rare variants or "ALL" if no exclusions',
            default => 'NMD_TRANSCRIPT,NON_SYNONYMOUS_CODING:NMD_TRANSCRIPT,NON_SYNONYMOUS_CODING,SPLICE_SITE:NMD_TRANSCRIPT,STOP_LOST:NON_SYNONYMOUS_CODING:NON_SYNONYMOUS_CODING,SPLICE_SITE:STOP_GAINED:STOP_GAINED,SPLICE_SITE'
        },
        select_phenotypes => {
            is => 'Text',
            is_many => 1,
            doc => "If specified, don't use all phenotypes from glm-model-file file, but instead only use these from a comma-delimited list -- this list's names must be an exact match to names specified in the glm-model-file",
            is_optional => 1
        },
        sample_list_file => {
            is => 'Text',
            doc => "Limit Samples in the Variant Matrix to Samples Within this File - Sample_Id should be the first column of a tab-delimited file, all other columns are ignored",
            is_optional => 1,
        },
        missing_value_markers => {
            is => 'Text',
            doc => 'Comma-delimited list of symbols that represent missing values such as "NA", ".", or "-"',
            is_optional => 1,
            default => 'NA,.,-999'
        },
        glm_model_file => {
            is => 'Text',
            doc => 'File outlining the type of model, response variable, covariants, etc. for the GLM analysis. (See DESCRIPTION).',
        },
        permutation_seed => {
            is => 'Text',
            doc => "If specified, uses the same permutation order for all analyses to keep correlations between traits",
            default => '123'
        },
        use_bsub => {
            is => "Boolean",
            doc => "If specified then Workflow is not used and instead thousands of jobs are bsub'd",
            default => 0
        },
        aggregate_only => {
            is => "Boolean",
            doc => "If specified then the output directory is assumed to contain completed results",
            default => 0
        },
        num_cores => {
            is => "Integer",
            doc => "The number of CPU cores to use in each gene/trait calculation",
            default => 4,
        },
    ],
};
sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run a burden analysis on germline (PhenotypeCorrelation) data"
}

sub help_synopsis {
    return <<EOS
Run a burden analysis on germline (PhenotypeCorrelation) data
EXAMPLE:    gmt germline burden-analysis --help
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
  return (
<<"EOS"

The clinical data file (glm-clinical-data-file) must follow these conventions:

=over 4

=item * Headers are required

=item * Each file must include at least 1 sample_id column and 1 attribute column, with the format being [sample_id  clinical_data_attribute_1 clinical_data_attribute_2 ...]

=item * The sample ID must match the sample ID listed in the VCF for relating the mutations of this sample. It uses an exact match, so all samples must have exactly the same nomenclature.

=back

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

Example:

/gscuser/qzhang/gstat/burdentest/readme    (readme file)
/gscuser/qzhang/gstat/burdentest/option_file_asms    (option file, this sets up the whole program)
/gscuser/qzhang/gstat/burdentest/burdentest.R  (main program)
/gscuser/qzhang/gstat/burdentest/rarelib.R (library )

/gscuser/qzhang/gstat/burdentest/jobs (job examples)
/gscmnt/sata424/info/medseq/Freimer-Boehnke/burdentest20120205/results (results)

EOS
    );
}


###############

sub execute {                               # replace with real execution logic.
    my $self = shift;

    $DB::single = 1;
    my $p_value_permutations = $self->p_value_permutations;
    my $seed = $self->permutation_seed;
    my $mutation_file = $self->mutation_file;
    my $phenotype_file = $self->glm_clinical_data_file;
    my $VEP_annotation_file = $self->VEP_annotation_file;
    my $output_directory = $self->output_directory;
    my $num_cores = $self->num_cores;

    my $project_name = $self->project_name;
    my $maf_cutoff = $self->maf_cutoff;
    my $permutations = $self->permutations;

    my %selected_phenotypes = map { $_ => 1 } $self->select_phenotypes;

    #parse the glm file to determine the clinical data types, names and covariates
    my $glm_model_file = $self->glm_model_file;
    my $glm_model_fh = Genome::Sys->open_file_for_reading($glm_model_file);
    my $glm_header = $glm_model_fh->getline;#This isn't used
    my %pheno_covar_type_hash;
    while(my $line = $glm_model_fh->getline) {
        next if ($line =~ m/^#/);
        chomp($line);
        my ($analysis_type,$clinical_data_trait_name,$variant_name,$covariates,$memo) = split(/\t/,$line);
        if($covariates eq 'NA') {
            $covariates = 'NONE';
        }
        if(!$self->select_phenotypes || exists($selected_phenotypes{$clinical_data_trait_name})) {
            $pheno_covar_type_hash{$clinical_data_trait_name} = "$covariates\t$analysis_type";
        }
    }
    $self->status_message("Phenotype/covariates selection: " . Dumper(\%pheno_covar_type_hash));

    #determine which trv_types to include
    #Note that the natural place to do this is now upstream in gmt vcf vcf-to-burden-matrix. This could still be useful for subsetting files created with that command though
    my $trv_types = $self->trv_types;
    my @trv_array = split(/:/, $trv_types);
    my $trv_types_to_use = q{"} . join(q{","}, @trv_array) . q{"};

    #set up missing value markers
    my $missing_value_markers = $self->missing_value_markers;
    my @missing_values = split(/,/,$missing_value_markers);
    my $R_missing_values = q{""};
    foreach my $marker (@missing_values) {
        $R_missing_values .= qq{,"$marker"};
    }

    #define subset of samples to use
    my %sample_name_hash;
    my $mutation_subset_file;
    if(defined($self->sample_list_file)) {
        my $sample_list_file = $self->sample_list_file;
        my $sample_list_inFh = Genome::Sys->open_file_for_reading($sample_list_file);
        while(my $sample_line = $sample_list_inFh->getline ) {
            chomp($sample_line);
            my ($sample_name) = split(/\t/, $sample_line);
            $sample_name_hash{$sample_name}++;
        }
        close($sample_list_inFh);

        $mutation_subset_file = "$output_directory/mutation_matrix.sample_subset.txt";
        my $fh_mutmat_out = Genome::Sys->open_file_for_overwriting($mutation_subset_file);
        unless ($fh_mutmat_out) {
            die "Failed to create new mutation matrix subset $mutation_subset_file!: $!";
        }
        my $mutmat_inFh = Genome::Sys->open_file_for_reading($mutation_file);
        my $mutmat_header = $mutmat_inFh->getline;
        my ($name, @samples) = split(/\t/, $mutmat_header);
        my %index_for_sample;
        @index_for_sample{@samples} = 0..$#samples;
        my @indexes_to_include = @index_for_sample{keys %sample_name_hash};

        print $fh_mutmat_out join("\t",@samples[@indexes_to_include]),"\n";

        while(my $line = $mutmat_inFh->getline ) {
            chomp($line);
            my ($variant_name, @values) = split(/\t/, $line);
            print $fh_mutmat_out join("\t",$variant_name,@values[@indexes_to_include]),"\n";
        }
        close($mutmat_inFh);
        close($fh_mutmat_out);
    }
    else {
        $mutation_subset_file = $mutation_file;
    }

    #get phenos and determine if trait is binary or not
    my $pheno_fh = Genome::Sys->open_file_for_reading($phenotype_file);
    my $pheno_header = $pheno_fh->getline;
    chomp($pheno_header);
    my @pheno_headers = split(/\t/, $pheno_header);
    my $subject_column_header = shift(@pheno_headers);
    $subject_column_header =~ s/ /\./g; #FIXME I don't think this really fixes the space issue. Wouldn't that be taken care of by R anyways? If not then Qunyuan could certainly add it.

    my $annot_fh = Genome::Sys->open_file_for_reading($VEP_annotation_file);
    my $annot_header = $annot_fh->getline;
    chomp($annot_header);
    my @header_fields = split /\t/, $annot_header;
    my %index_for_header;
    @index_for_header{@header_fields} = 0..$#header_fields;
    #NOTE that this is hardcoded for the VEP output.
    my $gene_name_in_header = 'Gene';

    my %gene_names;
    while (my $line = $annot_fh->getline) {
        chomp $line;
        my @fields = split(/\t/, $line);
        my $gene_name = $fields[$index_for_header{$gene_name_in_header}];
        next unless $gene_name; #this is to handle empty lines in the input file. Will no longer be needed once matrix generation is fixed.
        $gene_names{$gene_name} = 1;
    }

    #make .R option file
    my $R_option_file = "$output_directory/R_option_file.R";
    my $fh_R_option = Genome::Sys->open_file_for_overwriting($R_option_file);

    #-------------------------------------------------
    my $R_command_option = <<"_END_OF_R_";
### This is option file for burdentest.R ###

################### data files & key columns
missing.data=c($R_missing_values)

genotype.file="$mutation_subset_file"
gfile.delimiter="\\t"
gfile.vid="$project_name"   # variant id in genotype.file
gfile.sid="FIRST_ROW"              # subject id in genotype.file

phenotype.file="$phenotype_file"
pfile.delimiter="\\t"
pfile.sid="$subject_column_header"    # subject id in phenotype.file

anno.file="$VEP_annotation_file"
afile.delimiter="\\t"
afile.vid="Uploaded_variation"  #variant id in anno.file
gene.col="$gene_name_in_header"
vtype.col="Consequence"
vtype.use=c($trv_types_to_use)

out.dir="$output_directory"
if (!file.exists(out.dir)==T) dir.create(out.dir)

samplelist.dir="NONE" #right now this is not implemented in PhenotypeCorrelation so stub this in...

covariates="NONE"

########################### other options
maf.cutoff=$maf_cutoff

num.cores=$num_cores

_END_OF_R_
#covariates="$covariates"
    #-------------------------------------------------

    print $fh_R_option "$R_command_option\n";
    $fh_R_option->close;


    unless($self->aggregate_only) {
        #now create bsub commands
        #bsub -e err 'R --no-save < burdentest.R option_file_asms Q trigRES ABCA1 10000 123'



        #NOTE
        unless($self->use_bsub) {
            my %inputs;
            $inputs{option_file} = $R_option_file;
            $inputs{base_R_commands} = $base_R_commands;
            $inputs{permutations} = $self->permutations;
            $inputs{seed} = $self->permutation_seed;
            $inputs{p_value_permutations} = $self->p_value_permutations;

            my @outputs;

            my $job_number = 0;
            foreach my $phenotype (sort keys %pheno_covar_type_hash) {
                my ($covariates,$analysis_data_type) = split(/\t/,$pheno_covar_type_hash{$phenotype});
                foreach my $gene (keys %gene_names) {
                    $job_number++;
                    if ($gene eq '-' || $gene eq 'NA' ) {
                        next;
                    }
                    $inputs{"${job_number}_phenotype"} = $phenotype;
                    $inputs{"${job_number}_gene"} = $gene;
                    $inputs{"${job_number}_datatype"} = $analysis_data_type;
                    if($covariates && $covariates ne 'NONE') {
                        $inputs{"${job_number}_covariates"} = $covariates;
                    }
                }
            }
            my $workflow = Workflow::Model->create(
                name => 'Burden Test per gene per phenotype generation',
                input_properties => [
                keys %inputs
                ],
                output_properties => [
                map { 'test' . $_ . '_result' } 1..$job_number ,
                ],
            );
            for(my $job = 1; $job <= $job_number; $job++) {
                my $job_name = "burden test of " . $inputs{"${job}_gene"} . " for " . $inputs{"${job}_phenotype"};
                my $op = $workflow->add_operation(
                    name => $job_name,
                    operation_type => Workflow::OperationType::Command->get("Genome::Model::Tools::Germline::BurdenTest"),
                );
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"option_file",
                    right_operation=>$op,
                    right_property=>"option_file",
                );
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"permutations",
                    right_operation=>$op,
                    right_property=>"permutations",
                );
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"seed",
                    right_operation=>$op,
                    right_property=>"seed",
                );
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"p_value_permutations",
                    right_operation=>$op,
                    right_property=>"p_value_permutations",
                );

                #job specific inputs here
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"${job}_gene",
                    right_operation=>$op,
                    right_property=>"gene_name",
                );
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"${job}_datatype",
                    right_operation=>$op,
                    right_property=>"analysis_data_type",
                );
                $workflow->add_link(
                    left_operation=>$workflow->get_input_connector,
                    left_property=>"${job}_phenotype",
                    right_operation=>$op,
                    right_property=>"phenotype_name",
                );
                if(exists($inputs{"${job}_covariates"})) {
                    $workflow->add_link(
                        left_operation=>$workflow->get_input_connector,
                        left_property=>"${job}_covariates",
                        right_operation=>$op,
                        right_property=>"covariates",
                    );
                }

                #add output links
                $workflow->add_link(
                    left_operation=>$op,
                    left_property=>"result",
                    right_operation=>$workflow->get_output_connector,
                    right_property=>"test${job}_result",
                );
            }
            my @errors = $workflow->validate;
            $workflow->log_dir($self->output_directory);
            if (@errors) {
                $self->error_message(@errors);
                die "Errors validating workflow\n";
            }
            $self->status_message("Now launching burden test jobs");
            my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
            unless($result) {
                $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
                die $self->error_message("parallel burden test generation workflow did not return correctly.");
            }
        }
        else {
            #WORKFLOW does not scale well above so allow us to launch via bsub and hope for the best
            my $R_error_file = "$output_directory/R_error_file.err";
            my $bsub_base = "bsub -e $R_error_file 'Rscript $base_R_commands $R_option_file";
            foreach my $phenotype (sort keys %pheno_covar_type_hash) {
                my ($covariates,$analysis_data_type) = split(/\t/,$pheno_covar_type_hash{$phenotype});
                foreach my $gene (keys %gene_names) {
                    if ($gene eq '-' || $gene eq 'NA' ) {
                        next;
                    }
                    my $bsub_cmd;
                    if ($covariates && $covariates ne 'NONE') {
                        $bsub_cmd = "$bsub_base $analysis_data_type $phenotype $gene $permutations:$seed:$p_value_permutations $covariates\'";
                    }
                    else {
                        $bsub_cmd = "$bsub_base $analysis_data_type $phenotype $gene $permutations:$seed:$p_value_permutations\'";
                    }
                    $self->status_message("CMD: $bsub_cmd");
                    system($bsub_cmd);
                }
            }
        }
    }

    #NOTE When using bsub there's no way to know if all 20000 jobs have finished with the current implementation so rather than crash, just skip the summary
    if(!$self->use_bsub || $self->aggregate_only) {
        $self->status_message("Now aggregating burden test result files");
        my $header;
        my $num_fields;
        my @lines;  #since we won't know if we are NULL or not ahead of time, aggregate the header and lines together and then print
        my $header_source_file;
        my @null_pairs;
        foreach my $phenotype (sort keys %pheno_covar_type_hash) {
            my ($covariates,$analysis_data_type) = split(/\t/,$pheno_covar_type_hash{$phenotype});
            foreach my $gene (keys %gene_names) {
                #for each gene, find files, strip header and then aggregate result lines
                #if there is a null file then go ahead and enter a NULL line
                #what about error?
                my $file_stem = $self->output_directory . "/${phenotype}_${gene}";
                my $result_file = "$file_stem.burden.csv";
                my $null_file = "$file_stem.null";
                my $error_file = "$file_stem.error";
                if(-e $result_file) {
                    #we have a RESULT!
                    my $fh = Genome::Sys->open_file_for_reading($result_file);
                    #these should ALWAYS be two lines the way we are writing them, but let's assume the first line is the header and then other lines are results
                    my $header_line = $fh->getline;
                    chomp $header_line;
                    my @header_fields = split(",", $header_line);
                    if ($header) {
                        if (@header_fields != $num_fields) {
                            confess "Header in file $result_file does not match that in $header_source_file:\n"
                                ."$header_source_file:$header\n"
                                ."$result_file:$header_line"
                        }
                    } else {
                        $header = $header_line;
                        $header_source_file = $result_file;
                        $num_fields = @header_fields;
                    }
                    while(my $line = $fh->getline) {
                        chomp $line;
                        push @lines, $line if $line;    #hopefully this would catch empty lines.
                    }
                }
                elsif(-e $null_file) {
                    push @null_pairs, [$phenotype, $gene];
                }
                elsif(-e $error_file) {
                    #Workflow should have caught this, but lets just double check
                    $self->error_message("Error calculating burden test for $phenotype and gene $gene")
                }
                else {
                    #who knows what's up!?!
                    die $self->error_message("None of the expected output files for $phenotype and gene $gene produced");
                }
            }
        }

        my $ofh = Genome::Sys->open_file_for_overwriting($self->output_directory . "/burden_test_summary.csv");
        $ofh->print("$header\n");
        $ofh->print(join("\n",
            @lines, # non-NULL pheno/gene pairs
            (map { join(",", @$_, ("NA") x ($num_fields-2)) } @null_pairs) # NULL pheno/gene pairs
            ) . "\n"
        );
        $ofh->close;
    }

    return 1;
}

1;
