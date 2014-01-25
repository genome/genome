package Genome::Model::PhenotypeCorrelation::Command::CreateBurdenTestOptions;

use Genome;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::CreateBurdenTestOptions {
    is => "Command::V2",
    has_input => [
        mutation_file => {
            is => 'Text',
            doc => "Mutation Matrix",
        },
        glm_clinical_data_file => {
            is => 'Text',
            doc => "Phenotype File",
        },
        vep_annotation_file => {
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
        missing_value_markers => {
            is => 'Text',
            doc => 'A list of symbols that represent missing values such as "NA", ".", or "-"',
            is_optional => 1,
            is_many => 1,
            default => ["", "NA", ".", "-999"],
        },
        num_cores => {
            is => "Integer",
            doc => "The number of CPU cores to use in each gene/trait calculation",
            default => 4,
        },
    ],
    has_output => [
        output_file => {
            is => "FilePath",
            doc => "The output R file containing burden test options",
        }
    ],
};

sub execute {
    my $self = shift;
    my $output_path = $self->output_directory;
    my $output_file = "$output_path/R_option_file.R";
    my $missing_values = sprintf('"%s"', join('","', $self->missing_value_markers));

    my $clinical_data = Genome::Model::PhenotypeCorrelation::ClinicalData->from_file(
            $self->glm_clinical_data_file
            );

    my $subject_column_header = $clinical_data->subject_column_header;

    my $mutation_matrix = $self->mutation_file;
    my $project_name = $self->project_name;
    my $phenotype_file = $self->glm_clinical_data_file;
    my $vep_annotation_file = $self->vep_annotation_file;
    my $maf_cutoff = $self->maf_cutoff;
    my $num_cores = $self->num_cores;

    # Hard coded for VEP
    my $gene_name_in_header = 'Gene';

    my $out_fh = Genome::Sys->open_file_for_overwriting($output_file);
    $out_fh->write(<<EOF
### This is option file for burdentest.R ###

################### data files & key columns
missing.data=c($missing_values)

genotype.file="$mutation_matrix"
gfile.delimiter="\\t"
gfile.vid="$project_name"   # variant id in genotype.file
gfile.sid="FIRST_ROW"              # subject id in genotype.file

phenotype.file="$phenotype_file"
pfile.delimiter="\\t"
pfile.sid="$subject_column_header"    # subject id in phenotype.file

anno.file="$vep_annotation_file"
afile.delimiter="\\t"
afile.vid="Uploaded_variation"  #variant id in anno.file
gene.col="$gene_name_in_header"
vtype.col="Consequence"
vtype.use=c("ALL")

out.dir="$output_path"
if (!file.exists(out.dir)==T) dir.create(out.dir)

samplelist.dir="NONE" #right now this is not implemented in PhenotypeCorrelation so stub this in...

covariates="NONE"

########################### other options
maf.cutoff=$maf_cutoff

num.cores=$num_cores
EOF
);

    $out_fh->close;
    $self->output_file($output_file);
    return 1;
}

1;
