package Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile;

use strict;
use warnings;
use Genome;

class Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile {
    is => 'Command::V2',
    has_input => [
        tumor_rnaseq_builds => {
            is => 'Genome::Model::Build::RnaSeq',
            is_many => 1,
            is_input => 1,
        },
        normal_rnaseq_builds => {
            is => 'Genome::Model::Build::RnaSeq',
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
        },
        groups_list => {
            is => 'Text',
        }
    ],
};

sub execute {
    my $self = shift;

    my @tumor_builds = $self->tumor_rnaseq_builds;
    my @normal_builds = $self->normal_rnaseq_builds;

    my (@input_gene_count_files, @subjects, @groups);

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

    my $join_cmd = $self->_create_file_join_command($self->gene_counts_file, @input_gene_count_files);

    my $rv = Genome::Sys->shellcmd(cmd => $join_cmd);

    unless ($rv) {
        die $self->error_message("Join command unsuccessful");
    }

    $self->subjects_list(join ",", @subjects);
    $self->groups_list(join ",", @groups);

    return 1;
}

sub _retrieve_build_information {
    my $self = shift;
    my $group_name = shift;
    my @builds = @_;

    my (@input_gene_count_files, @groups, @subjects);

    for my $build (@builds) {
        push @input_gene_count_files, $build->data_directory . "/results/digital_expression_result/gene-counts.tsv";
        push @subjects, $build->subject->source;
        push @groups, $group_name
    }

    return { input_gene_count_files => \@input_gene_count_files, subjects => \@subjects, groups => \@groups };
}

sub _create_file_join_command {
    my $self = shift;
    my $output_file = shift;
    my @files = @_;

    if (scalar(@files) < 2) {
        die $self->error_message("You need at least two files to join");
    }

    my $cmd = "join " . $files[0] . " " . $files[1];

    $cmd = join (" \| join - ", $cmd, @files[2..$#files]) . " > " . $output_file;

    return $cmd;
}

1;
