package Genome::Model::RnaSeq::Command::AnnotateChimerascan;

use strict;
use warnings;
use Genome;
use File::Spec;

class Genome::Model::RnaSeq::Command::AnnotateChimerascan {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Genome::Model::Build::RnaSeq',
            is_output => 1,
        },
        cancer_annotation_db_id => {
            is => 'Text',
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::RnaSeq',
            id_by => 'build_id',
        },
        cancer_annotation_db => {
            is => 'Genome::Db::Tgi::CancerAnnotation',
            id_by => 'cancer_annotation_db_id',
            doc => 'db of cancer annotation (see \'genome db list\' for latest version of desired database)',
        },
        annotated_bedpe_file => {
            is => 'Path',
            is_output => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    my $cancer_annotation_db = $self->cancer_annotation_db;
    my $infile = $self->_infile;
    my @gene_name_columns = $self->_gene_name_columns;
    my $gene_groups = $self->_gene_groups;

    my $cmd = Genome::Model::ClinSeq::Command::AnnotateGenesByCategory->create(
        infile => $infile,
        cancer_annotation_db => $cancer_annotation_db,
        gene_name_columns => \@gene_name_columns,
        gene_groups => $gene_groups,
    );

    my $rv = $cmd->execute;
    unless($rv){
        $self->error_message("Failed to annotate Chimerascan output: $@");
        return;
    }

    $self->annotated_bedpe_file($cmd->category_outfile);
    
    return 1;
}

sub _infile {
    my $self = shift;
    my $build = $self->build;
    return File::Spec->join($build->data_directory, 'fusions', 'filtered_chimeras.bedpe');
}

sub _gene_name_columns {
    my $self = shift;
    return qw/ 3P 5P /;
}

sub _gene_groups {
    return 'Fusions';
}

1;
