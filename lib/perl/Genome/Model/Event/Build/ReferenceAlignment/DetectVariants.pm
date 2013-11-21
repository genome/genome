package Genome::Model::Event::Build::ReferenceAlignment::DetectVariants;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::ReferenceAlignment::DetectVariants{
    is => ['Genome::Model::Event'],
};

sub bsub_rusage {
    my $self = shift;

    return "-R 'select[tmp>4000] rusage[tmp=4000]' -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER}";
}

sub execute{
    my $self = shift;

    $self->status_message("Executing detect variants step");
    my $build = $self->build;

    my $processing_profile = $build->processing_profile;

    my %params;
    $params{snv_detection_strategy} = $processing_profile->snv_detection_strategy if $processing_profile->snv_detection_strategy;
    $params{indel_detection_strategy} = $processing_profile->indel_detection_strategy if $processing_profile->indel_detection_strategy;
    $params{sv_detection_strategy} = $processing_profile->sv_detection_strategy if $processing_profile->sv_detection_strategy;
    $params{cnv_detection_strategy} = $processing_profile->cnv_detection_strategy if $processing_profile->cnv_detection_strategy;

    my $bam = $build->whole_rmdup_bam_file;
    $params{aligned_reads_input} = $bam;

    my $reference_build = $build->model->reference_sequence_build;
    my $reference_fasta = $reference_build->full_consensus_path('fa');
    unless(-e $reference_fasta){
        die $self->error_message("fasta file for reference build doesn't exist!");
    }
    $params{reference_build_id} = $reference_build->id;

    my $output_dir = $build->variants_directory;
    $params{output_directory} = $output_dir;

    my $aligned_reads_sample = $build->subject->name;
    $params{aligned_reads_sample} = $aligned_reads_sample;

    my $command = Genome::Model::Tools::DetectVariants2::Dispatcher->create(%params);
    unless ($command){
        die $self->error_message("Couldn't create detect variants dispatcher from params:\n".Data::Dumper::Dumper \%params);
    }
    my $rv = $command->execute;
    my $err = $@;
    unless ($rv){
        die $self->error_message("Failed to execute detect variants dispatcher(err:$@) with params:\n".Data::Dumper::Dumper \%params);
    }
    else {
        my @results = $command->results;
        my $test_name = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || '';
        push @results, map { Genome::Model::Tools::DetectVariants2::Result::Vcf->get(input_id => $_->id, test_name => $test_name); } @results;
        for my $result (@results) {
            $result->add_user(user => $build, label => 'uses');
        }
    }

    $self->status_message("detect variants command completed successfully");

    return 1;
}

1;
