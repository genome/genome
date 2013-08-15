package Genome::Model::SomaticValidation::Command::SubmissionSummary;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::SomaticValidation::Command::SubmissionSummary {
    is => 'Command::V2',
    doc => "List a summary of the merged alignment BAMs for the provided builds and a file suitable for submitting the bam list.",
    has => [
        sample_mapping_file => {
            is=> 'String',
            doc => 'this command will generate this list of bam file names & samples',
        }, 
        bam_list_file => {
            is=> 'String',
            doc => 'this command will generate this list of bam file names, newline separated suitable for gxfer to submit bams',
        },
        md5_list_file => {
            is=> 'String',
            doc => 'this will generate a list of bam md5 file names, newline separated',
        },
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'List of models to use'
        },
    ],
};


sub help_detail {
    return "List the path of the merged alignment BAMs for the provided builds/flow cells, and generate a sample mapping.";
}


sub execute {
    my $self = shift;

    my @builds = map{$_->last_succeeded_build} $self->models;

    my $samp_map = IO::File->new(">".$self->sample_mapping_file);
    unless ($samp_map) {
        $self->error_message("Failed to open sample mapping file for writing ". $self->sample_mapping_file);
        return;
    }
    my $bam_list = IO::File->new(">".$self->bam_list_file);
    unless ($bam_list) {
        $self->error_message("Failed to open bam list file for writing ". $self->bam_list_file);
        return;
    }

    my $md5_list = IO::File->new(">".$self->md5_list_file);
    unless ($md5_list) {
        $self->error_message("Failed to open md5 list file for writing ". $self->md5_list_file);
        return;
    }

    for my $build (@builds) {
        my $roi_name = $build->model->region_of_interest_set_name ? $build->model->region_of_interest_set_name : 'N/A';
        my $refbuild_name = $build->model->reference_sequence_build->name ? $build->model->reference_sequence_build->name : 'N/A';

        print $samp_map join ("\t", $build->tumor_sample->name, $refbuild_name, $roi_name, $build->tumor_bam, basename($build->tumor_bam));
        print $samp_map "\n";
        if($build->normal_sample) {
            print $samp_map join ("\t", $build->normal_sample->name, $refbuild_name, $roi_name, $build->normal_bam, basename($build->normal_bam));
            print $samp_map "\n";
        }
    }

    my @bams = map {$_->tumor_bam} @builds;
    push @bams, grep { defined $_ } map {$_->normal_bam} @builds;
    print $bam_list join ("\n", @bams),"\n";
    print $md5_list join ("\n", map {$_.".md5"} @bams), "\n";
    for my $md5 (map {$_.".md5"} @bams) {
        print STDERR "$md5 doesn't exist\n" if ! -e $md5;
    }

    $samp_map->close;
    $bam_list->close;
    $md5_list->close;

    return 1;
}

1;

