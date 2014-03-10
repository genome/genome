package Genome::Model::RnaSeq::Command::MakeGeneCountsFile;

use strict;
use warnings;
use Genome;

class Genome::Model::RnaSeq::Command::MakeGeneCountsFile {
    is => 'Command::V2',
    has_input => [
        tumor_models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            is_input => 1,
        },
        normal_models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            is_input => 1,
        },
    ],
    has_output => [
        counts_file => {
            is => 'File',
        },
        subjects => {
            is => 'Text',
            is_optional => 1,
        },
        groups => {
            is => 'Text',
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;

    my (@tumor_builds, @normal_builds);

    for my $model ($self->tumor_models) {
        if ( defined($model->last_succeeded_build) ) {
                push @tumor_builds, $model->last_succeeded_build;
        }
    }
    for my $model ($self->normal_models) {
        if ( defined($model->last_succeeded_build) ) {
                push @normal_builds, $model->last_succeeded_build;
        }
    }

    _validate_all_inputs_have_same_annotation(@normal_builds, @tumor_builds);

    my (@input_counts_files, @header, @subjects, @groups);

    my $tumor_builds_information = $self->_retrieve_build_information("tumor", @tumor_builds);
    my $normal_builds_information = $self->_retrieve_build_information("normal", @normal_builds);

    push @input_counts_files,
        @{$tumor_builds_information->{input_counts_files}},
        @{$normal_builds_information->{input_counts_files}} ;
    push @subjects,
        @{$tumor_builds_information->{subjects}},
        @{$normal_builds_information->{subjects}} ;
    push @groups,
        @{$tumor_builds_information->{groups}},
        @{$normal_builds_information->{groups}} ;
    push @header,
        @{$tumor_builds_information->{headers}},
        @{$normal_builds_information->{headers}};

    $self->_check_gene_column_identical(@input_counts_files);

    $self->_write_header_to_file($self->counts_file, @header);

    my $join_cmd = $self->_create_file_join_command($self->counts_file, @input_counts_files);

    my $rv = Genome::Sys->shellcmd(cmd => $join_cmd);

    unless ($rv) {
        die $self->error_message("Join command unsuccessful");
    }

    $self->subjects(join ",", @subjects);
    $self->groups(join ",", @groups);

    return 1;
}

sub _validate_all_inputs_have_same_annotation {
    my $self = shift;
    my @builds = @_;

    my $annot_build;

    for my $build (@builds) {
        unless (defined $build->annotation_build) {
            die $self->error_message("Build ".$build->id." does not have annotation build defined");
        }
        unless (defined $annot_build) {
            $annot_build = $build->annotation_build->id;
        }
        unless ($build->annotation_build->id eq $annot_build) {
            die $self->error_message("Builds must have the same annotation input");
        }
    }

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

    my (@input_counts_files, @groups, @subjects, @header);

    for my $build (@builds) {
        my $subject = $build->subject->source->common_name;
        push @input_counts_files, $build->data_directory . "/results/digital_expression_result/gene-counts.tsv";
        push @subjects, $subject;
        push @groups, $group_name;
        push @header, "$subject-$group_name";
    }

    return { input_counts_files => \@input_counts_files, 
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
