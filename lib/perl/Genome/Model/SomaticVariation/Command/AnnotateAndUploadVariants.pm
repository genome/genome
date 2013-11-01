package Genome::Model::SomaticVariation::Command::AnnotateAndUploadVariants;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Genome::Model::Tools::DetectVariants2::Utilities qw(
    final_result_for_variant_type
);

class Genome::Model::SomaticVariation::Command::AnnotateAndUploadVariants{
    is => 'Genome::Command::Base',
    has =>[
        build_id => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticVariation model',
        },
        annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run   during the annotation step',
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4 ],
            is_optional => 1,
            is_input => 1,
        },
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
            id_by => 'build_id',
        },
        regulatory_annotations => {
            is => 'Genome::FeatureList',
            is_input => 1,
            is_many => 1,
            is_optional => 1,
        },
        get_regulome_db => {
            is => 'Boolean',
            default => 0,
            is_optional => 1,
            is_input => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
        lsf_resource => {
            default => "-M 20000000 -R 'select[mem>20000] rusage[mem=20000]'",
        },
        joinx_version => {
            is => "Text",
            doc => "Version of joinx to use",
            is_optional => 0,
            default => 1.6,
        },
        lsf_resource => { 
            is => 'Text',
            default => "-R 'rusage[mem=16000]' -M 16000000",
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

    my ($tier1_snvs, $tier2_snvs, $tier3_snvs, $tier4_snvs, $tier1_indels, $tier2_indels, $tier3_indels, $tier4_indels);
    
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
        $tier3_snvs = $build->data_set_path("effects/snvs.hq.novel.tier3", $version, "bed");
        unless(-e $tier3_snvs) {
            die $self->error_message("No tier 3 snvs file for build!");
        }
        $files{'snvs.hq.tier3'} = $tier3_snvs;
        $tier4_snvs = $build->data_set_path("effects/snvs.hq.novel.tier4", $version, "bed");
        unless(-e $tier4_snvs) {
            die $self->error_message("No tier 4 snvs file for build!");
        }
        $files{'snvs.hq.tier4'} = $tier4_snvs;
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
        $tier3_indels = $build->data_set_path("effects/indels.hq.novel.tier3",$version,'bed');
        unless(-e $tier3_indels) {
            die $self->error_message("No tier 3 indels file for build!");
        }
        $files{'indels.hq.tier3'} = $tier3_indels;
        $tier4_indels = $build->data_set_path("effects/indels.hq.novel.tier4",$version,'bed');
        unless(-e $tier4_indels) {
            die $self->error_message("No tier 4 indels file for build!");
        }
        $files{'indels.hq.tier4'} = $tier4_indels;
    }

    #annotate variants
    my $annotator_version = $self->annotator_version;
    unless($annotator_version){
        die $self->error_message("No variant annotator version for build!");
    }

    my $annotation_build_id = $self->build->annotation_build->id;

    my %annotation_params = (
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
                } 
                else {
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
                
                #Make a version of the annotation file without a header
                `mv $final_output_file $final_output_file.header`;
                `grep -v "^chromosome_name" $final_output_file.header > $final_output_file`;
            }
        }
        else{
            $self->status_message("No variants present for $key, skipping annotation and upload");
            File::Copy::copy($variant_file, $annotated_file);
            File::Copy::copy($variant_file, "$annotated_file.top");
            File::Copy::copy($variant_file, $uploaded_file);
        }
    }

    #Annotate vcfs with dbsnp ids
    $self->status_message("Adding dbsnp ids to vcf");
    
    if ($build->previously_discovered_variations_build) {
        my $annotation_vcf = $build->previously_discovered_variations_build->snvs_vcf;
        if ($annotation_vcf && -e $annotation_vcf) {
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
    
    #Annotate SVs. Need to handle mouse vs human and human 36 vs 37
    $self->status_message("Executing sv annotation");

    my $species_name = $build->subject->species_name;
    
    if ($build->sv_detection_strategy) {
        my $sv_sr  = final_result_for_variant_type([$build->results], 'sv');
        if ($sv_sr) {
            my $sv_dir = $sv_sr->output_dir;
            if ($sv_dir and -d $sv_dir) {
                my $sv_file = $sv_dir . '/svs.merge.file.somatic';
                if ($sv_file and -s $sv_file) {
                    if ($species_name and $species_name =~/^(human|mouse)$/) {
                        my $sv_annot_out_file = $build->data_directory . '/effects/svs.hq.annotated';
                        my $fusion_out_file   = $build->data_directory . '/effects/svs.hq.fusion_transcripts.out';
                        my @annotator_list    = qw(Transcripts FusionTranscripts);

                        my $cancer_gene_list = join("/",Genome::Sys->dbpath(join("/","cancer-gene-list",$species_name), 1),"Cancer_genes.csv");

                        my %params = (
                            input_file  => $sv_file,
                            output_file => $sv_annot_out_file,
                            fusion_transcripts_fusion_output_file => $fusion_out_file,
                            annotation_build_id => $annotation_build_id,
                            annotator_list      => \@annotator_list,
                            transcripts_print_flanking_genes => 1,
                            transcripts_cancer_gene_list     => $cancer_gene_list,
                            chrA_column       => 1,
                            bpA_column        => 2,
                            chrB_column       => 4,
                            bpB_column        => 5,
                            event_type_column => 7,
                            orient_column     => 8,
                            score_column      => 12,  
                        );
                        
                        my $ref_seq_build = $build->reference_sequence_build;
                        if ($species_name eq 'human') {
                            my $annot_file = $self->_get_human_sv_annot_file();

                            if ($annot_file) {
                                push @annotator_list, 'Dbsnp', 'Segdup', 'Dbvar';
                                %params = (
                                    %params,
                                    annotator_list => \@annotator_list,
                                    dbsnp_annotation_file  => $annot_file->{dbsnp},
                                    segdup_annotation_file => $annot_file->{segdup},
                                    dbvar_annotation_file  => $annot_file->{dbvar},
                                );
                            }
                            else {
                                $self->warning_message('Human annotation for dbsnp, dbvar and segdup is only set for build 36 and 37');
                            }
                        }
                        elsif ($species_name eq 'mouse') {
                            unless ($ref_seq_build->id == 107494762) {#UCSC-mouse build37, mouse ref seq used in production
                                $self->warning_message('For now sv annotation only run on UCSC-mouse build37');
                                %params = ();
                            }
                        }

                        if (%params) {
                            my $sv_annot = Genome::Model::Tools::Annotate::Sv::Combine->create(%params);
                            my $rv = $sv_annot->execute;    
                            $self->warning_message("sv annotation did not finish ok") unless $rv == 1;
                        }
                    }
                    else {
                        $self->warning_message('The build is not for human or mouse. Skip sv annotation');
                    }
                }
                else {
                    $self->warning_message("combine sv somatic file: $sv_file is not valid");
                }
            }
            else {
                $self->warning_message("sv dir: $sv_dir is not valid");
            }
        }
        else {
            $self->warning_message('No sv software result for this build');
        }
    }

    if ($self->regulatory_annotations) {
        my $count = 0;
        my $final_name = $build->data_directory."/effects/snvs.hq.annotated.bed";
        my $header = "#chr\tstart\tstop\talleles\tcount1\tcount2";
        my @added_columns;
        my $input_file = $build->data_directory."/variants/snvs.hq.bed";
        for my $annotation ($self->regulatory_annotations) {
            if ($annotation->is_1_based) {
                $self->error_message("Regulatory annotation sets must be in true-BED format");
                die $self->error_message;
            }
            $count++;
            push @added_columns, $count+6;
            my $base_name = $annotation->name;
            my $bed_file = $annotation->file_path;
            my $output_file = $build->data_directory."/effects/snvs.hq.".$count.".bed";
            my $rt = Genome::Model::Tools::BedTools::Map->execute(
                input_file_a => $input_file,
                input_file_b => $bed_file,
                output_file => $output_file,
                operation => "distinct",
                column => 4,
                use_version => "2.16.2",
                null => '-',
            );
            unless ($rt) {
                $self->error_message("Failed to annotate with bedtools map");
                return;
            }
            $header .= "\t$base_name";
            $input_file = $output_file;
        }
        
        my $cmd = "echo \"$header\" > $final_name; cat ".$build->data_directory."/effects/snvs.hq.$count.bed >> $final_name";
        $self->status_message("Running command $cmd");
        `$cmd`;

        my @in_files;
        for my $tier ((1,2,3,4)) {
            if(my $file = $build->data_set_path("effects/snvs.hq.tier$tier", $annotation_output_version, "annotated.top.header")){
                if (-s $file) {
                    push @in_files, $file;
                }
            }
        }
        for my $in_file (@in_files) {
            my $append = Genome::Model::Tools::Annotate::AppendColumns->execute(
                additional_columns_file => $final_name,
                input_variants => $in_file,
                output_file => $in_file.".regulatory",
                columns_to_append => \@added_columns,
                chrom_column => 1,
                start_column => 2,
                stop_column => 3,
            );
            unless ($append) {
                $self->error_message("Append columns failed for ".$in_file);
                return;
            }
        }
    }
    if ($self->get_regulome_db) {
        my $rdb_file = $build->data_set_path("effects/snvs.hq.regulomedb", 1, "full");
        my $rdb_rv = Genome::Model::Tools::RegulomeDb::GetAnnotationsForVariants->execute(
            variant_list => $build->data_set_path("variants/snvs.hq", 2, "bed"),
            output_file => $rdb_file,
            format => "full",
        );
    }

    #upload variants

    ##upload metrics
    $self->status_message("Uploading metrics");
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

sub _get_db_version {
    my $self = shift;
    my $ref_seq_build = $self->build->reference_sequence_build;
    my $version = $ref_seq_build->version;
    my $type;
    if ($version) {
        if ($version =~ /36/) {
            $type = "build36";
        }
        elsif ($version =~ /37/) {
            $type = "build37";
        }
    } 
    if ($ref_seq_build->id == 109104543) {
        $type = "build36";
    }   
    return $type;
}

#only run on human build 36 and 37 for now
sub _get_human_sv_annot_file {
    my $self = shift;
    my $type = $self->_get_db_version();
    my $dbsnp_version;
    my $species = "human";

    if ($type eq "build36") {
        $dbsnp_version = "130";
    }
    elsif ($type eq "build37") {
        $dbsnp_version = "132";
    }

    return unless $type and $dbsnp_version;
    $self->status_message("Getting annotation files for species $species of version $type, dbsnp version $dbsnp_version");
    my $segdup = Genome::Db->get(source_name => 'ucsc', database_name => $species, external_version => $type);
    my $dbsnp = Genome::Db->get(source_name => 'genome-db-dbsnp', database_name => "$species/$type", external_version => $dbsnp_version);
    my $dbvar = Genome::Db->get(source_name => 'dbvar', database_name => $species, external_version => $type);
    return {
        segdup => $segdup->data_directory . '/segdup.tsv',
        dbsnp  => $dbsnp->data_directory . '/dbsnp.csv',
        dbvar  => $dbvar->data_directory . '/dbvar.tsv',
        };
}

1;
