package Genome::Model::SomaticVariation::Command::AnnotateAndUploadVariants;

use strict;
use warnings;
use Genome;
use Data::Dumper;

class Genome::Model::SomaticVariation::Command::AnnotateAndUploadVariants{
    is => 'Genome::Command::Base',
    has =>[
        build_id => {
            is => 'Integer',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticVariation model',
        },
        annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run   during the annotation step',
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3],
            is_optional => 1,
            is_input => 1,
        },
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
            id_by => 'build_id',
        }
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
    ],
};

sub execute{
    my $self = shift;
    my $build = $self->build;
    unless ($build){
        die $self->error_message("no build provided!");
    }

    $self->status_message("Executing Annotate and Upload step");

    my $version = 2;
    #my $version = GMT:BED:CONVERT::version();  TODO, something like this instead of hardcoding
    
    my %files;

    my ($tier1_snvs, $tier2_snvs, $tier1_indels, $tier2_indels);
    
    if ($build->snv_detection_strategy){
        $tier1_snvs = $build->data_set_path("effects/snvs.hq.novel.tier1",$version,'bed');
        unless(-e $tier1_snvs){
            die $self->error_message("No tier 1 snvs file for build!");
        }
        $files{'snvs.hq.tier1'} = $tier1_snvs;

        $tier2_snvs = $build->data_set_path("effects/snvs.hq.novel.tier2",$version,'bed');
        unless(-e $tier2_snvs){
            die $self->error_message("No tier 2 snvs file for build!");
        }
        $files{'snvs.hq.tier2'} = $tier2_snvs;
    }

    if ($build->indel_detection_strategy){
        $tier1_indels = $build->data_set_path("effects/indels.hq.novel.tier1",$version,'bed');
        unless(-e $tier1_indels){
            die $self->error_message("No tier 1 indels file for build!");
        }
        $files{'indels.hq.tier1'} = $tier1_indels;

        $tier2_indels = $build->data_set_path("effects/indels.hq.novel.tier2",$version,'bed');
        unless(-e $tier2_indels){
            die $self->error_message("No tier 2 indels file for build!");
        }
        $files{'indels.hq.tier2'} = $tier2_indels;
    }

    #annotate variants
    my $annotator_version = $self->annotator_version;
    unless($annotator_version){
        die $self->error_message("No variant annotator version for build!");
    }

    my $annotation_build_id = $self->build->annotation_build->id;

    my %annotation_params = (
        no_headers => 1,
        use_version => $annotator_version,
        build_id => $annotation_build_id,
    );

    my $annotation_output_version=1;  #TODO, do we even have an annotation file format?  need to figure out how to resolve this w/ data_set_path
    my $upload_output_version = 1; #TODO, same issue here as with annotation output version

    for my $key (keys %files){
        my $variant_file = $files{$key};
        unless ($variant_file){
            $self->status_message("No detection strategy for $key. Skipping annotation and upload");
        }
        unless (-e $variant_file){
            die $self->error_message("File expected for annotating $key at $variant_file does not exist!  Failing");
        }

        my $annotated_file = $build->data_set_path("effects/$key",$annotation_output_version,'annotated');
        my $uploaded_file = $build->data_set_path("effects/$key", $upload_output_version, "uploaded");

        $annotation_params{variant_bed_file} = $variant_file;
        
        if (-s $variant_file){

            # Run annotation twice, once for the "none" annotation filter and once for the "top" annotation filter
            for my $annotation_filter ( qw(none top) ) {
                $annotation_params{annotation_filter} = $annotation_filter;

                # "None" will the be the default named output file, since annotated.none sounds a bit weird.
                if ($annotation_filter eq "none") {
                    $annotation_params{output_file} = $annotated_file;
                } else {
                    $annotation_params{output_file} = "$annotated_file.$annotation_filter";
                }

                $self->status_message("Creating TranscriptVariants object with parameters: " . Dumper \%annotation_params);
                my $annotation_command = Genome::Model::Tools::Annotate::TranscriptVariants->create(%annotation_params);
                unless ($annotation_command){
                    die $self->error_message("Failed to create annotate command for $key. Params:\n".Data::Dumper::Dumper(\%annotation_params));
                }
                my $annotate_rv = $annotation_command->execute;
                my $annotate_err = $@;
                unless($annotate_rv){
                    die $self->error_message("Failed to execute annotate command for $key(err: $annotate_err) from params:\n" . Data::Dumper::Dumper(\%annotation_params));
                }
                unless(-s $annotation_params{output_file}){
                    die $self->error_message("No output from annotate command for $key. Params:\n" . Data::Dumper::Dumper(\%annotation_params));
                }
            }

            my %upload_params = (
                build_id => $self->build_id, 
            );

            $upload_params{variant_file} = $variant_file;
            $upload_params{annotation_file} = $annotated_file;
            $upload_params{output_file} = $uploaded_file;

            my $upload_command = Genome::Model::Tools::Somatic::UploadVariants->create(%upload_params);
            unless ($upload_command){
                die $self->error_message("Failed to create upload command for $key. Params:\n". Data::Dumper::Dumper(\%upload_params));
            }
            #my $upload_rv = $upload_command->execute;
            my $upload_rv =  1;  #TODO turn this on
            my $upload_err = $@;
            unless ($upload_rv){
                die $self->error_message("Failed to execute upload command for $key (err: $upload_err). Params:\n".Dumper(\%upload_params));
            }
            unless(-s $upload_params{output_file} or 1){
                die $self->error_message("No output from upload command for $key. Params:\n" . Data::Dumper::Dumper(\%upload_params));
            }
        }else{
            $self->status_message("No variants present for $key, skipping annotation and upload");
            File::Copy::copy($variant_file, $annotated_file);
            File::Copy::copy($variant_file, "$annotated_file.top");
            File::Copy::copy($variant_file, $uploaded_file);
        }
    }
    
    #upload variants

    ##upload metrics
    eval{
        my ($tier3_snvs, $tier4_snvs, $tier3_indels, $tier4_indels, $svs);
        $tier3_snvs = $build->data_set_path("effects/snvs.hq.novel.tier3",$version,'bed');
        $tier4_snvs = $build->data_set_path("effects/snvs.hq.novel.tier4",$version,'bed');
        $tier3_indels = $build->data_set_path("effects/indels.hq.novel.tier3",$version,'bed');
        $tier4_indels = $build->data_set_path("effects/indels.hq.novel.tier4",$version,'bed');

        $svs = $build->data_set_path("variants/svs.hq",$version,'bed');

        my ($number_of_tier1_snvs, $number_of_tier2_snvs, $number_of_tier3_snvs, $number_of_tier4_snvs, $number_of_tier1_indels, $number_of_tier2_indels, $number_of_tier3_indels, $number_of_tier4_indels, $number_of_svs);
        if($build->snv_detection_strategy) {
            chomp($number_of_tier1_snvs = `wc -l $tier1_snvs | cut -f 1 -d' '`);
            chomp($number_of_tier2_snvs = `wc -l $tier2_snvs | cut -f 1 -d' '`);
            chomp($number_of_tier3_snvs = `wc -l $tier3_snvs | cut -f 1 -d' '`);
            chomp($number_of_tier4_snvs = `wc -l $tier4_snvs | cut -f 1 -d' '`);
            $build->set_metric("tier1_snv_count", $number_of_tier1_snvs);
            $build->set_metric("tier2_snv_count", $number_of_tier2_snvs);
            $build->set_metric("tier3_snv_count", $number_of_tier3_snvs);
            $build->set_metric("tier4_snv_count", $number_of_tier4_snvs);
        }


        if($build->indel_detection_strategy) {
            chomp($number_of_tier1_indels = `wc -l $tier1_indels | cut -f 1 -d' '`);
            chomp($number_of_tier2_indels = `wc -l $tier2_indels | cut -f 1 -d' '`);
            chomp($number_of_tier3_indels = `wc -l $tier3_indels | cut -f 1 -d' '`);
            chomp($number_of_tier4_indels = `wc -l $tier4_indels | cut -f 1 -d' '`);
            $build->set_metric("tier1_indel_count", $number_of_tier1_indels);
            $build->set_metric("tier2_indel_count", $number_of_tier2_indels);
            $build->set_metric("tier3_indel_count", $number_of_tier3_indels);
            $build->set_metric("tier4_indel_count", $number_of_tier4_indels);
        }

        if($build->sv_detection_strategy) {
            chomp($number_of_svs = `wc -l $svs | cut -f 1 -d' '`); 
            $build->set_metric("sv_count", $number_of_svs);
        }
    }; 
    $self->status_message("Metric setting problem: $@") if $@;

     
    return 1;
}

1;
