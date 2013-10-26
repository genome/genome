package Genome::Model::PhenotypeCorrelation::Command::ParallelClinicalCorrelation;

use Sort::Naturally qw/nsort/;
use File::Basename qw/basename/;
use File::Temp;
use Genome;
use Workflow::Simple;
use POSIX qw/WIFEXITED/;
use Carp qw/confess/;
use File::Path qw/mkpath/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::ParallelClinicalCorrelation {
    is => ["Genome::Command::Base"],
    doc => "Run the clinical correlation tool from the music suite.",
    has_input => [
        variant_matrix => {
            is => "File",
            doc => "Variant matrix to analyze",
        },
        work_directory => {
            is => "File",
            doc => "Where to store the intermediate output files. This must be a network location visible from all blades.",
        },
        cleanup_intermediate_files => {
            is => "Boolean",
            doc => "If true, intermediate files are removed (unset for debugging)",
            default_value => 1,
        },
        output_prefix => {
            is => "File",
            doc => "Where to write the final merged result",
            is_output => 1,
        },
        sample_list_file => {
            is => "String",
            doc => 'File containing samples names, 1 per line, for input into MuSiC',
        },
        clinical_data_file => {
            is => "String",
            doc => 'File containing clinical data',
        },
        glm_model_file => {
            is => "String",
            doc => 'File containing the model specification for this analysis',
        },
        max_cols_per_file => {
            is => "Number",
            doc => "Max number of columns per submatrix",
        },
    ],
    has_output => [
        categorical_output_file => {
            is => "Text",
            doc => "The path to the glm output file",
            is_optional => 1,
        },
        glm_output_file => {
            is => "Text",
            doc => "The path to the glm output file",
            is_optional => 1,
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
This command splits the input variant matrix into files containing at most
\$max_cols_per_file columns each and processes them with music's clinical
correlation command in parallel.
EOS
}

sub execute {
    my $self = shift;

    my $max_cols = $self->max_cols_per_file;
    my $glm_model_file = $self->glm_model_file;
    my $clin_file = $self->clinical_data_file;
    my $samples_file = $self->sample_list_file;
    my $variant_matrix = $self->variant_matrix;
    my $cleanup = $self->cleanup_intermediate_files || 0;

    $self->status_message("Cleanup is set to $cleanup");

    my $work_dir_obj = File::Temp->newdir(
        'ParallelClinicalCorrelation-intermediate-XXXXX',
        DIR => $self->work_directory,
        CLEANUP => $cleanup,
    );
    my $work_dir = $work_dir_obj->dirname;

    my $log_dir = Workflow::Model->parent_workflow_log_dir;
    unless($log_dir) {
        $log_dir = $work_dir . "/logs";
        unless(-d $log_dir) {
            unless(mkpath($log_dir)) {
                die $self->error_message("Unable to create log directories for correlations");
            }
        }
    }

    $self->work_directory($work_dir);

    my $glm_output_file = $self->output_prefix . ".glm.tsv";
    my $categorical_output_file = $self->output_prefix . ".categorical.tsv";
    my $submatrix_dir = "$work_dir/submatrices";
    my $sub_results_dir = "$work_dir/sub_results";

    my $categorical_data_file = "$work_dir/categorical_clinical_data.txt";
    $self->_make_categorical_data_file($categorical_data_file);

    confess "max_cols_per_file must be > 0" unless $max_cols > 0;
    confess "Output directory $work_dir does not exist" unless -d $work_dir;
    for my $d ($submatrix_dir, $sub_results_dir) {
        mkdir($d);
        confess "Failed to create directory $d" unless -d $d;
    }

    my $split_cmd = Genome::Model::PhenotypeCorrelation::Command::SplitVariantMatrix->create(
        input_file => $variant_matrix,
        output_prefix => "$submatrix_dir/submatrix",
        max_cols_per_file => $max_cols
    );
    my $submatrices = $split_cmd->execute;
    confess "Failed to split variant matrix $variant_matrix!" unless @$submatrices;

    $self->status_message("Split $variant_matrix into " . scalar(@$submatrices) . " submatrices, creating workflow.");

    if ($ENV{WF_USE_FLOW}) {
        $self->status_message("Sleeping 90 seconds for NFS cache.");
        sleep(90);
    }

    my %params = (
        sample_list_file => $samples_file,
        output_directory => $sub_results_dir,
        glm_model_file => $glm_model_file,
        clinical_data_file => $clin_file,
        variant_matrix => $submatrices
    );
    $params{categorical_clinical_data_file} = $categorical_data_file if -s $categorical_data_file;

    my $op = Workflow::Operation->create(
        name => "Parallel clinical correlation",
        operation_type => Workflow::OperationType::Command->get(
            'Genome::Model::PhenotypeCorrelation::Command::ClinicalCorrelationWrapper',
            )
    );
    $op->parallel_by('variant_matrix');
    $op->log_dir($log_dir);

    $self->status_message("Executing workflow...");
    my $output = Workflow::Simple::run_workflow_lsf($op, %params);
    unless (defined $output) {
        $self->error_message("Workflow execution failed!");
        my @error;
        for (@Workflow::Simple::ERROR) {
            push @error, $_->error;
        }
        $self->error_message(join("\n", @error));
        confess $self->error_message;
    }
    $self->status_message("Workflow completed, merging intermediate results...");

    my @glm_results = nsort glob("$sub_results_dir/*.glm.tsv");
    my @categorical_results = nsort glob("$sub_results_dir/*.categorical.tsv");

    if (@glm_results) {
        $self->status_message("Merging " .scalar(@glm_results). " glm results");
        my $merged_glm_fh = Genome::Sys->open_file_for_writing($glm_output_file);
        $self->_merge_results($merged_glm_fh, @glm_results);
        $self->glm_output_file($glm_output_file);
    }
    else {
        $self->status_message("No glm results produced");
    }
    if (@categorical_results) {
        $self->status_message("Merging " .scalar(@categorical_results). " categorical results");
        my $tempfile = "$categorical_output_file.tmp";
        my $merged_categorical_fh = Genome::Sys->open_file_for_writing($tempfile);
        $self->_merge_results($merged_categorical_fh, @categorical_results);
        $merged_categorical_fh->close();
        $self->_recalculate_categorical_stats($tempfile, $categorical_output_file);
        unlink($tempfile);
        $self->categorical_output_file($categorical_output_file);
    }
    else {
        $self->status_message("No categorical results produced");
    }


    $self->status_message("Parallel Clinical Correlation completed successfully.");

    return 1;
};

sub _make_categorical_data_file {
    my ($self, $output_file) = @_;
    my $glm_model_file = $self->glm_model_file;
    my $clin_file = $self->clinical_data_file;
    my $glm_config = Genome::Model::PhenotypeCorrelation::GlmConfig->from_file($glm_model_file);

    my @categorical_attrs = map {$_->{attr_name}} $glm_config->categorical_attributes;

    my $out_fh = Genome::Sys->open_file_for_overwriting($output_file);
    my $clinical_data = Genome::Model::PhenotypeCorrelation::ClinicalData->from_file($clin_file);
    $clinical_data->to_filehandle($out_fh, attribute_names => \@categorical_attrs);
}

sub _merge_results {
    my ($self, $ofh, @files) = @_;

    # input/output file handles
    my @ifhs = map {Genome::Sys->open_file_for_reading($_)} @files;
    $self->status_message("Opened " . scalar(@ifhs) . " intermediate files for merging.");

    my $header = $ifhs[0]->getline();
    for my $i (1..$#ifhs) {
        my $ifh = $ifhs[$i];
        my $xheader = $ifh->getline();
        if ($xheader ne $header) {
            confess "Multiple header formats found while merging clinical "
                ."correlation results:\n\t"
                .join("\n\t", $header, $xheader);
        }
    }

    $ofh->write($header);
    for my $i (0..$#ifhs) {
        $self->status_message("Processing file $i of " . scalar(@ifhs));
        my $ifh = $ifhs[$i];
        while (my $line = $ifh->getline()) {
            $ofh->write($line);
        }
    }
}

# The categorical test does p-value correction (FDR and Bonferroni) which depends on
# having all of the p-values. Since we're working in parallel, these values are
# computed incorrectly for each chunk. Once the file has been merged, we can load it
# into R to recalculate the correct values.
sub _recalculate_categorical_stats {
    my ($self, $infile, $outfile) = @_;

    my $R_code = <<EOR;
tt <- read.table("$infile", sep="\t", header=TRUE);
tt\$FDR = p.adjust(tt\$Pval, "fdr");

# Match formatting in GMT/Music/ClinicalCorrelation.pm.R
tt[,"Statistic"] = sapply(tt[,"Statistic"], sprintf, fmt="%.4E");
tt[,"Pval"] = sapply(tt[,"Pval"], sprintf, fmt="%.4E");
tt[,"FDR"] = sapply(tt[,"FDR"], sprintf, fmt="%.2E");

# Match ordering in GMT/Music/ClinicalCorrelation.pm.R
tt=tt[order(tt[,"Gene"]),];
tt=tt[order(tt[,"Pval"]),];

colnames(tt)=c("Gene","ClinParam","Method","NumCases","Statistic","Pval","FDR");

write.table(tt,file="$outfile",quote=FALSE,row.names=FALSE,sep="\t");

EOR
    my $R_script = "$outfile.tmp.R";
    my $fh = Genome::Sys->open_file_for_overwriting($R_script);
    $fh->write($R_code);
    $fh->close();
    my $R_cmd = "R --slave --args < $R_script";
    my $rv = system $R_cmd;
    unlink($R_script);
    WIFEXITED($rv) or confess "Couldn't run: $R_cmd ($rv)";
}

1;
