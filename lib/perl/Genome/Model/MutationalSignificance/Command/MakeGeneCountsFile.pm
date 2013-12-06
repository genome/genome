package Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile;

use strict;
use warnings;
use Genome;

class Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile {
    is => 'Command::V2',
    has_input => [
        clinseq_models => {
            is => 'Genome::Model::ClinSeq',
            is_many => 1,
            is_input => 1,
        },
    ],
    has_output => [
        gene_counts_file => {
            is => 'File',
        },
        subjects_list => {
            is => 'Text',
            is_optional => 1,
        },
        groups_list => {
            is => 'Text',
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;

    my @clinseq_models = $self->clinseq_models;

    my (@tumor_builds, @normal_builds);

    for my $model (@clinseq_models) {
        if ( defined($model->tumor_rnaseq_model) && defined($model->tumor_rnaseq_model->last_succeeded_build) ) {
                push @tumor_builds, $model->tumor_rnaseq_model->last_succeeded_build;
        }
        if ( defined($model->normal_rnaseq_model) && defined($model->normal_rnaseq_model->last_succeeded_build) ) {
                push @normal_builds, $model->normal_rnaseq_model->last_succeeded_build;
        }
    }

    my (@input_gene_count_files, @header, @subjects, @groups);

    my $tumor_builds_information = $self->_retrieve_build_information("tumor", @tumor_builds);
    my $normal_builds_information = $self->_retrieve_build_information("normal", @normal_builds);

    push @input_gene_count_files,
        @{$tumor_builds_information->{input_gene_count_files}},
        @{$normal_builds_information->{input_gene_count_files}} ;
    push @subjects,
        @{$tumor_builds_information->{subjects}},
        @{$normal_builds_information->{subjects}} ;
    push @groups,
        @{$tumor_builds_information->{groups}},
        @{$normal_builds_information->{groups}} ;
    push @header,
        @{$tumor_builds_information->{headers}},
        @{$normal_builds_information->{headers}};

    $self->_check_gene_column_identical(@input_gene_count_files);

    $self->_write_header_to_file($self->gene_counts_file, @header);

    my $join_cmd = $self->_create_file_join_command($self->gene_counts_file, @input_gene_count_files);

    my $rv = Genome::Sys->shellcmd(cmd => $join_cmd);

    unless ($rv) {
        die $self->error_message("Join command unsuccessful");
    }

    $self->subjects_list(join ",", @subjects);
    $self->groups_list(join ",", @groups);

    return 1;
}

sub _write_header_to_file {
    my $self = shift;
    my $output_file = shift;
    my @headers = @_;

    my $out = Genome::Sys->open_file_for_writing($output_file);
    $out->print(join("\t", "GENE", @headers)."\n");
    $out->close;

    return 1;
}

sub _retrieve_build_information {
    my $self = shift;
    my $group_name = shift;
    my @builds = @_;

    my (@input_gene_count_files, @groups, @subjects, @header);

    for my $build (@builds) {
        my $subject = $build->subject->source->common_name;
        push @input_gene_count_files, $build->data_directory . "/results/digital_expression_result/gene-counts.tsv";
        push @subjects, $subject;
        push @groups, $group_name;
        push @header, "$subject-$group_name";
    }

    return { input_gene_count_files => \@input_gene_count_files, 
             subjects => \@subjects, 
             groups => \@groups, 
             headers => \@header };
}

sub _check_gene_column_identical {
    my $self = shift;
    my @files = @_;

    if (scalar(@files) < 2) {
        die $self->error_message("You need at least two files to join");
    }

    my $cut_cmd = "cut -f1 " . $files[0];
    my $master_column = `$cut_cmd`;

    for my $file (@files) {
         $cut_cmd = "cut -f1 " . $file;
         my $comparison_column = `$cut_cmd`;

         unless ($comparison_column eq $master_column) {
            die $self->error_message("Gene column of $file is not the same");
         }
    }

    return 1;
}

sub _create_file_join_command {
    my $self = shift;
    my $output_file = shift;
    my @files = @_;

    if (scalar(@files) < 2) {
        die $self->error_message("You need at least two files to join");
    }

    my $cmd = "join -t '\t' " . $files[0] . " " . $files[1];

    $cmd = join (" \| join -t '\t' - ", $cmd, @files[2..$#files]) . " >> " . $output_file;

    return $cmd;
}

1;
