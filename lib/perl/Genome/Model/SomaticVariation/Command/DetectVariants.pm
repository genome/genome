package Genome::Model::SomaticVariation::Command::DetectVariants;

use strict;
use warnings;
use Genome;

class Genome::Model::SomaticVariation::Command::DetectVariants{
    is => 'Genome::Command::Base',
    has =>[
        build_id => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticVariation model',
        },
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
            id_by => 'build_id',
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
};

sub execute{
    my $self = shift;

    $self->debug_message("Executing detect variants step");
    my $build = $self->build;
    unless ($build){
        die $self->error_message("no build provided!");
    }

    my %params;
    $params{snv_detection_strategy} = $build->snv_detection_strategy if $build->snv_detection_strategy;
    $params{indel_detection_strategy} = $build->indel_detection_strategy if $build->indel_detection_strategy;
    $params{sv_detection_strategy} = $build->sv_detection_strategy if $build->sv_detection_strategy;
    $params{cnv_detection_strategy} = $build->cnv_detection_strategy if $build->cnv_detection_strategy;
    
    my $tumor_bam = $build->tumor_bam;
    unless (-e $tumor_bam){
        die $self->error_message("No tumor bam found for somatic model");
    }
    $params{aligned_reads_input} = $tumor_bam;

    my $reference_build = $build->reference_sequence_build;
    my $reference_fasta = $reference_build->full_consensus_path('fa');
    unless(-e $reference_fasta){
        die $self->error_message("fasta file for reference build doesn't exist!");
    }
    $params{reference_build_id} = $reference_build->id;
    
    my $normal_bam = $build->normal_bam;
    unless (-e $normal_bam){
        die $self->error_message("No normal bam found for somatic model");
    }
    $params{control_aligned_reads_input} = $normal_bam;

    my $output_dir = $build->data_directory."/variants";
    $params{output_directory} = $output_dir;

    my $tmodel = $build->tumor_model;
    my $nmodel = $build->normal_model;
    my $aligned_reads_sample = $tmodel->subject->name;
    my $control_aligned_reads_sample = $nmodel->subject->name;
    $params{aligned_reads_sample} = $aligned_reads_sample;
    $params{control_aligned_reads_sample} = $control_aligned_reads_sample;
    
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
        my @results = $command->results, $command->lq_results;
        push @results, map { Genome::Model::Tools::DetectVariants2::Result::Vcf->get(input_id => $_->id ); } @results;
        for my $result (@results) {
            $result->add_user(user => $build, label => 'uses');
        }
    }

    $self->debug_message("detect variants command completed successfully");

    my $version = 2;
    #my $version = GMT:BED:CONVERT::version();  TODO, something like this instead of hardcoding

    if ($build->snv_detection_strategy){
        my $result = $build->data_set_path("variants/snvs.hq",$version,'bed'); 
        unless (-e $result){
            my $unexpected_format_output = $command->_snv_hq_output_file;
            unless (-e $unexpected_format_output){
                die $self->error_message("Expected hq detected snvs file $result, but it does not exist!");
            }
            symlink($unexpected_format_output, $result);
        }
        my $lq_result = $build->data_set_path("variants/snvs.lq",$version,'bed');
        unless (-e $lq_result){
            #these are not set on the detect variants command, so I'm hardcoding them.
            my $unexpected_filename_output = $build->data_directory."/variants/snvs.lq.bed";
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected lq detected snvs $unexpected_filename_output, but it does not exist");
            }
            symlink($unexpected_filename_output, $lq_result);
        }

    }

    if ($build->indel_detection_strategy){
        my $result = $build->data_set_path("variants/indels.hq",$version,'bed'); 
        unless (-e $result){
            my $unexpected_filename_output = $command->_indel_hq_output_file;
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected hq detected indels file $result, but it does not exist!");
            }
            symlink($unexpected_filename_output, $result);
        }
        my $lq_result = $build->data_set_path("variants/indels.lq",$version,'bed');
        unless (-e $lq_result){
            #these are not set on the detect variants command, so I'm hardcoding them.
            my $unexpected_filename_output = $build->data_directory."/variants/indels.lq.bed";
            if (-e $unexpected_filename_output){
                symlink($unexpected_filename_output, $lq_result);
            } else {
                $self->debug_message("No lq indel file found. Creating an empty file at $lq_result");
                system("touch $lq_result");
            }
        }
    }

    if ($build->sv_detection_strategy){
        my $result = $build->data_set_path("variants/svs.hq",$version,'bed'); 
        unless (-e $result){
            my $unexpected_filename_output = $command->_sv_hq_output_file;
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected hq detected snvs file $result, but it does not exist!");
            }
            symlink($unexpected_filename_output, $result);
        }
    }

    if ($build->cnv_detection_strategy){
        my $result = $build->data_set_path("variants/cnvs.hq",$version,'bed'); 
        unless (-e $result){
            my $unexpected_filename_output = $command->_cnv_hq_output_file;
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected hq detected snvs file $result, but it does not exist!");
            }
            symlink($unexpected_filename_output, $result);
        }
    }

    $self->_create_tcga_vcfs;

    $self->debug_message("detect variants step completed");

    return 1;
}

sub _create_tcga_vcfs {
    my $self = shift;
    my $build = $self->build;

    for my $type ("snv", "indel") {
        # Create a TCGA compliant vcf
        my $original_vcf = $build->data_directory . "/variants/" . $type . "s.detailed.vcf.gz";
        my $tcga_vcf = $build->data_directory . "/variants/" . $type . "s_tcga.tar.gz";
        my $strategy_accessor = $type . "_detection_strategy";
        if ($build->$strategy_accessor and -s $original_vcf) {
            my $tcga_command = Genome::Model::Tools::Vcf::TcgaSanitize->create(
                input_file => $original_vcf,
                output_file => $tcga_vcf,
                package_for_tcga => 1,
            );
            unless ($tcga_command->execute) {
                die $self->error_message("Failed to TcgaSanitize the final snvs vcf file");
            }
        }
    }

    return 1;
}

1;

