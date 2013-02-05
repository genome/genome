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
        joinx_version => {
            is => "Text",
            doc => "Version of joinx to use",
            is_optional => 0,
            default => 1.6,
        },
        lsf_resource => { 
            is => 'Text',
            default => "-R 'rusage[mem=8000]' -M 8000000",
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
    my %vcf_files;

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
        my $snvs_vcf = $build->data_directory."/variants/snvs.vcf.gz";
        unless(-e $snvs_vcf) {
            die $self->error_message("No snvs vcf");
        }
        $vcf_files{'snvs.hq'} = $snvs_vcf;
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
                my $final_output_file;
                if ($annotation_filter eq "none") {
                    $final_output_file = $annotated_file;
                    $annotation_params{output_file} = $final_output_file.".non-dedup";
                } else {
                    $final_output_file = "$annotated_file.$annotation_filter";
                    $annotation_params{output_file} = "$final_output_file.non-dedup";
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
                #Deduplicate the annotation file
                my $in = Genome::Sys->open_file_for_reading($annotation_params{output_file});
                my $output = Genome::Sys->open_file_for_writing($final_output_file);
                my %seen;
                while(my $line = <$in>) {
                    chomp $line;
                    if ($seen{$line}) {
                        next;
                    }
                    $seen{$line} = 1;
                    print $output $line."\n";
                }
                $in->close;
                $output->close;
                my $rm_cmd = "rm -f ".$annotation_params{output_file};
                `$rm_cmd`;
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

    #Annotate vcfs with dbsnp ids
    if ($build->previously_discovered_variations_build) {
        my $annotation_vcf = $build->previously_discovered_variations_build->snvs_vcf;
        if (-e $annotation_vcf) {
            for my $key (keys(%vcf_files)) {
                my $variant_file = $vcf_files{$key};
                my $output_file = $variant_file.".annotated.vcf.gz";
                $output_file =~ s/vcf.gz.//;
                my $info_string = $build->processing_profile->vcf_annotate_dbsnp_info_field_string eq "NO_INFO" ? "" : $build->processing_profile->vcf_annotate_dbsnp_info_field_string;
                my $info = $info_string eq "" ? 0 : 1;
                my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->execute(
                    input_file=> $variant_file,
                    annotation_file=>$annotation_vcf,
                    output_file=>$output_file,
                    use_bgzip=>1,
                    info_fields=>$info_string,
                    info => $info,
                    use_version => $self->joinx_version,
                ) || die "Failed to execute Joinx Vcf annotation using db: $annotation_vcf";
                $self->status_message("Successfully annotated VCF with information from $annotation_vcf");
            }
            foreach my $annotation_prefix ("snvs.hq.tier1", "snvs.hq.tier2") {
                my $annotation_top_file = $build->data_set_path("effects/$annotation_prefix", $annotation_output_version, "annotated.top");
                next unless (-e $annotation_top_file);
                my $vcf_file = $vcf_files{"snvs.hq"};
                $vcf_file = $vcf_file.".annotated.vcf.gz";
                $vcf_file =~ s/vcf.gz.//;
                $self->status_message("Adding rsid column for $annotation_prefix from vcf file ".$vcf_file);
                my $rv = Genome::Model::Tools::Annotate::AddRsid->execute(
                    anno_file => $annotation_top_file,
                    vcf_file => $vcf_file,
                    output_file => $build->data_set_path("effects/$annotation_prefix.rsid", $annotation_output_version, "annotated.top"),
                );
            }
        }
        else {
            $self->warning_message("No snvs vcf available for previously_discovered_variations_build");
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
